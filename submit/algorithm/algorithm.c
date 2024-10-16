#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "algorithm.h"
#include "LKHInterface.h"
#pragma GCC optimize ("O2")

#define DEBUG_LEVEL 0
// #define DEBUG_LEVEL 1
// #define DEBUG_LEVEL 2

#define INF INT_MAX/200

double MyGetTime(){ // 返回实际时间：秒
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}
double fabs(double number) {
    uint64_t mask = 0x7FFFFFFFFFFFFFFF; // 清除符号位的掩码
    uint64_t bits = *(uint64_t*)&number; // 将 double 转换为 uint64_t
    bits &= mask; // 清除符号位
    return *(double*)&bits; // 将结果转换回 double
}
double sqrt(double n) {
    if (n < 0) {
        // 对于负数没有实数平方根
        return -1;
    }
    if (n == 0) {
        return 0;
    }

    double x = n;  // 初始猜测值
    double tolerance = 1e-7;  // 误差容限
    double diff;

    do {
        double next_x = 0.5 * (x + n / x);
        diff = fabs(next_x - x);  // 计算与前一次的差值
        x = next_x;
    } while (diff > tolerance);

    return x;
}
double ceil(double x) {
    int int_part = (int)x; // 将 x 转换为整数部分

    // 如果 x 是正数并且不等于其整数部分，则需要向上取整
    if (x > int_part) {
        return int_part + 1.0;
    }

    // 如果 x 是负数或整数，则直接返回整数部分
    return int_part;
}

// 设置需要调整的 LKH 的参数
void loadUserChangedParam(int matDimension, LKHParameters *p, double scheduleStartTime){
    if(matDimension <= 1002){ // 10,50,100,1000 
        // 默认值
    }
    else if(matDimension <= 2002){ // 2000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        p->POPMUSIC_SampleSize = 20;
        p->POPMUSIC_Solutions = 50;
    }
    else if(matDimension <= 5002){ // 5000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        p->POPMUSIC_SampleSize = 20;
        p->POPMUSIC_Solutions = 20;
    }
    else{ // 10000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        p->POPMUSIC_SampleSize = 20;
        p->POPMUSIC_Solutions = 10;
    }
    p->Runs = 1;
    p->TraceLevel = 1;
    p->TimeLimit = DBL_MAX; // 由总时间、一定时间跨度内的改进值共同控制退出即可
    p->TotalTimeLimit = 400; // 最大允许运行时间
    p->ScheduleScoreInSecond = 20;
    p->MoveType = 3;
    p->TimeSpan = 2;
}

// LKH 算法
int32_t LKH(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    /* 生成邻接矩阵 */
    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    for (uint32_t i = 0; i < len; ++i){
        adjMat[i] = (int *)malloc(sizeof(int) * len);
    }
   
    for (uint32_t i = 0; i < len - 2; ++i) {//为了节点能和id号对应，还是将它们从0开始对应行列
        for (uint32_t j = 0; j < len - 2; ++j) {
            if(i == j)
                adjMat[i][j] = 0;
            else {
                HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                adjMat[i][j] = SeekTimeCalculate(&start, &end);
            }
        }
    }

     /*len-2行len-2列是虚拟节点，len-1行len-1列是磁头节点，为虚拟节点和磁头节点设置代价 */
    for (uint32_t i = 0; i < len - 2; ++i) {
        adjMat[len-2][i] = INF;  // 虚拟节点到其他节点的代价为极大值，则不可能达到其他点
        adjMat[i][len-2] = 0;  // 其他节点到虚拟节点的代价为0
        
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = SeekTimeCalculate(&start, &end);  // 磁头节点到其他节点的代价为对应寻址时间
        adjMat[i][len-1] = INF;  // 其他节点到磁头节点的代价为极大值，则不可能返回该点
    }
    adjMat[len-2][len-1] = 0; // 虚拟节点到磁头节点的代价为0，到其他结点代价无穷，相当于固定了磁头节点为第二个节点
    adjMat[len-2][len-2] = 0; // 虚拟节点到自身的代价为0
    adjMat[len-1][len-1] = 0; // 磁头节点到自身的代价为0
    adjMat[len-1][len-2] = INF;  // 磁头节点到虚拟节点的代价为极大值

    /* 打印邻接矩阵 */
#if DEBUG_LEVEL >= 2
    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            if (adjMat[i][j] == INF) {
                printf("INF\t");  // 打印无穷大（不可达）的情况
            } else {
                printf("%u\t ", adjMat[i][j]);  // 打印代价
            }
        }
        printf("\n");
    }
#endif

    /* 基于最近邻贪心构造初始解， 注意巡回路径从1开始编号*/
    int *intialTour = (int *)malloc(len * sizeof(int));
    int* visited = (int*) calloc(len + 1, sizeof(int)); // 标记已经到达过的城市，数组初始化为0
    // 从虚拟结点出发, 再到磁头结点
    intialTour[0] = len - 1, intialTour[1] = len; 
    visited[len-1] = 1, visited[len] = 1;
    int current = len; // 磁头
    for (int i = 2; i < len; ++i){
        int min_dist = INF;
        int node = -1;
        // 找到最近的未访问城市
        for (int j = 1; j <= len; ++j){
            if(!visited[j] && adjMat[current-1][j-1] < min_dist){
                min_dist = adjMat[current-1][j-1];
                node = j;
            }
        }
        // 更新当前城市
        visited[node] = 1;
        intialTour[i] = node;
        current = node;
    }

    /* 设置固定边, 结点从1开始编号 */
    int fixLen = 1;
    int *fixEdge = (int *)malloc(2 * fixLen * sizeof(int));
    // 虚拟结点到磁头结点的边固定
    fixEdge[0] = len - 1;
    fixEdge[1] = len;

    /* 调用 LKH 求解 */
    /* 确定LKH输入结构体 */
    LKHInput *lkhInput = (LKHInput *)malloc(sizeof(LKHInput));
    lkhInput->adjMat = adjMat;
    lkhInput->matDimension = len;
    lkhInput->intialTour = intialTour;
    lkhInput->fixEdge = fixEdge;
    lkhInput->fixEdgeLen = fixLen;
    lkhInput->scheduleStartTime = scheduleStartTime;
    lkhInput->lkhParameters = (LKHParameters *)malloc(sizeof(LKHParameters));
    loadDefaultParam(lkhInput->lkhParameters);
    loadUserChangedParam(lkhInput->matDimension, lkhInput->lkhParameters, scheduleStartTime);
    /* 确定LKH输出结构体 */
    LKHOutput *lkhOutput = (LKHOutput *)malloc(sizeof(LKHOutput));
    lkhOutput->tourCost = 0;
    lkhOutput->tourResult = (int *)malloc(len * sizeof(int));

    ret = solveTSP(lkhInput, lkhOutput);
   
    // 处理解
    int i, j;
    for (i = 0; i < len; ++i){
        if(lkhOutput->tourResult[i] == len) // i指向磁头结点
            break;
    }
    for (j = 0; j < len - 2; ++j){ // 排序结果赋值到输出
        if(i + 1 > len - 1)
            i = -1;
        output->sequence[j] = lkhOutput->tourResult[++i];
    }

    /* 打印求解结果：输出output->sequence的内容 */
#if DEBUG_LEVEL >= 1
    printf("LKH result\n");
    printf("io len:%d\n",len-2);
    printf("Output sequence:\n");
    for (uint32_t i = 0; i < len - 2; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放内存
    for (int i = 0; i < len; i++) {
        free(adjMat[i]);  // 逐行释放
    }
    free(adjMat);
    free(intialTour);
    free(visited);
    free(fixEdge);
    free(lkhInput->lkhParameters);
    free(lkhInput);
    free(lkhOutput->tourResult);
    free(lkhOutput);

    /* 调用公共函数示例：调用电机寻址、带体磨损、电机磨损函数 */
    // HeadInfo start = {input->ioVec.ioArray[0].wrap, input->ioVec.ioArray[0].endLpos, HEAD_RW};
    // HeadInfo end = {input->ioVec.ioArray[1].wrap, input->ioVec.ioArray[1].endLpos, HEAD_RW};
    // int32_t seekT = 0;
    // int32_t beltW = 0;
    // int32_t motorW = 0;
    //    for (uint32_t i = 0; i < 10000; i++) {
    //        seekT = SeekTimeCalculate(&start, &end);
    //        beltW = BeltWearTimes(&start, &end, NULL);
    //        motorW = MotorWearTimes(&start, &end);
    //    }

    // /* 调用公共函数示例：调用IO读写时间函数 */
    // uint32_t rwT = ReadTimeCalculate(abs(input->ioVec.ioArray[0].endLpos - input->ioVec.ioArray[0].startLpos));

    return ret;
}

// Sort 算法的自定义比较函数：首先按 wrap 排序，相同 wrap 的按 startLpos 排序
int Sort_compare(const void *a, const void *b) {
    IOUint *ioA = (IOUint *)a;
    IOUint *ioB = (IOUint *)b;

    if (ioA->wrap != ioB->wrap) {
        return (ioA->wrap > ioB->wrap) - (ioA->wrap < ioB->wrap);  // 按 wrap 升序
    } else {
        return (ioA->startLpos > ioB->startLpos) - (ioA->startLpos < ioB->startLpos);  // 按 startLpos 升序
    }
}

int32_t Sort(const InputParam *input, OutputParam *output){
    // 创建 ioArray 的副本，用于排序
    IOUint *sort_IOArray = (IOUint *)malloc(input->ioVec.len * sizeof(IOUint));
    
    // 将 input 的 ioArray 拷贝到 sortedIOArray 中
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        sort_IOArray[i] = input->ioVec.ioArray[i];
    }

    // 使用 qsort 对副本进行排序
    qsort(sort_IOArray, input->ioVec.len, sizeof(IOUint), Sort_compare);
    // 将排序后的 IOUint 的 id 存入 output->sequence 中
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        output->sequence[i] = sort_IOArray[i].id;
    }
#if DEBUG_LEVEL >= 1
    printf("Sort算法得到的io序列:");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放副本内存
    free(sort_IOArray);
    return 0;
}

// Scan 算法的自定义比较函数：按 startLpos 排序
int Scan_compare(const void *a, const void *b) {
    IOUint *ioA = (IOUint *)a;
    IOUint *ioB = (IOUint *)b;

    if (ioA->startLpos < ioB->startLpos) {
        return -1; // ioA 在 ioB 之前
    } else if (ioA->startLpos > ioB->startLpos) {
        return 1;  // ioA 在 ioB 之后
    } else {
        return 0;  // 两者相等
    }
}

int32_t Scan(const InputParam *input, OutputParam *output){    
    // 创建 ioArray 的副本，用于排序
    IOUint *scan_IOArray = (IOUint *)malloc(input->ioVec.len * sizeof(IOUint));

    // 将 input 的 ioArray 拷贝到 scan_IOArray 中
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        scan_IOArray[i] = input->ioVec.ioArray[i];
    }

    // 使用 qsort 对副本进行排序
    qsort(scan_IOArray, input->ioVec.len, sizeof(IOUint), Scan_compare);

    uint32_t seqIndex = 0;
    // 扫描方向从 BOT -> EOT (从小到大)
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        if (scan_IOArray[i].wrap % 2 == 0) {  // 假设 wrap 偶数表示从 BOT 向 EOT
            output->sequence[seqIndex++] = scan_IOArray[i].id;
        }
    }

    // 到达 EOT 后，反向扫描（从大到小）
    for (int i = input->ioVec.len - 1; i >= 0; --i) {
        if (scan_IOArray[i].wrap % 2 == 1) {  // 假设 wrap 奇数表示从 EOT 向 BOT
            output->sequence[seqIndex++] = scan_IOArray[i].id;
        }
    }
#if DEBUG_LEVEL >= 1
    printf("Scan算法得到的io序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
         printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放副本内存
    free(scan_IOArray);
    return 0;
}

// 最近邻贪心构造
int32_t NearestNeighbor(const InputParam *input, OutputParam *output){
    /* 生成邻接矩阵 */
    uint32_t len = output->len + 1; //增加了一个磁头起始位置节点
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    for (uint32_t i = 0; i < len; ++i){
        adjMat[i] = (int *)malloc(sizeof(int) * len);
    }
   
    for (uint32_t i = 0; i < len - 1; ++i) {//为了节点能和id号对应，还是将它们从0开始对应行列
        for (uint32_t j = 0; j < len - 1; ++j) {
            if(i == j)
                adjMat[i][j] = 0;
            else {
                HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                adjMat[i][j] = SeekTimeCalculate(&start, &end);
            }
        }
    }

     /*len-1行len-1列是磁头节点，为磁头节点设置代价 */
    for (uint32_t i = 0; i < len - 1; ++i) {  
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = SeekTimeCalculate(&start, &end);  // 磁头节点到其他节点的代价为对应寻址时间
        adjMat[i][len-1] = INF;  // 其他节点到磁头节点的代价为极大值，则不可能返回该点
    }
    adjMat[len-1][len-1] = 0; // 磁头节点到自身的代价为0

    /* 打印邻接矩阵 */
#if DEBUG_LEVEL >= 2
    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            if (adjMat[i][j] == INF) {
                printf("INF\t");  // 打印无穷大（不可达）的情况
            } else {
                printf("%u\t ", adjMat[i][j]);  // 打印代价
            }
        }
        printf("\n");
    }
#endif

    /* 贪心构造, 注意巡回路径从1开始编号 */
    int *tour = (int *)malloc(len * sizeof(int));
    int* visited = (int*) calloc(len + 1, sizeof(int)); // 标记已经到达过的城市，数组初始化为0
    // 从磁头出发
    tour[0] = len; 
    visited[len] = 1;
    int current = len;
    for (int i = 1; i < len; ++i){
        int min_dist = INF;
        int node = -1;
        // 找到最近的未访问城市
        for (int j = 1; j <= len; ++j){
            if(!visited[j] && adjMat[current-1][j-1] < min_dist){
                min_dist = adjMat[current-1][j-1];
                node = j;
            }
        }
        // 更新当前城市
        visited[node] = 1;
        tour[i] = node;
        current = node;
    }
    
    // 赋值到输出数据结构
    for (int i = 0; i < len - 1; ++i){ 
        output->sequence[i] = tour[i+1];
    }
#if DEBUG_LEVEL >= 1
    printf("最近邻算法得到的io序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
         printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif    

    // 释放内存
    for (int i = 0; i < len; i++) {
        free(adjMat[i]);  // 逐行释放
    }
    free(adjMat);
    free(tour);
    free(visited);

    return 0;
}

/**
 * @brief  算法接口
 * @param  input            输入参数
 * @param  output           输出参数
 * @return int32_t          返回成功或者失败，RETURN_OK 或 RETURN_ERROR
 */
int32_t IOScheduleAlgorithm(const InputParam *input, OutputParam *output)
{    
    // 使用 LKH 算法
    return LKH(input, output);
    // 使用最近邻贪心算法 
    // return NearestNeighbor(input, output);
    // 使用 Sort 算法
    // return Sort(input, output);
    // 使用 Scan 算法
    // return Scan(input, output);
}

/**
 * @brief  算法运行的主入口
 * @param  input            输入参数
 * @param  output           输出参数
 * @return uint32_t          返回成功或者失败，RETURN_OK 或 RETURN_ERROR
 */
int32_t AlgorithmRun(const InputParam *input, OutputParam *output)
{
    int32_t ret;

    ret = IOScheduleAlgorithm(input, output);

    return ret;
}

//################################### LKHInterface.c begin ###################################
#include "LKH.h"
#include "Genetic.h"
#include "Heap.h"
#include "LKHInterface.h"

/*
 * LKH内核 接口，欲调用LHK，按以下三个步骤走：
 * (1)设定 LKH内核 的相关参数值：函数 ReadParameters()。一般情况下无需手动设置，采用默认值即可;
 * (2)设定问题类型、规模等，并进行样例数据的输入：函数 ReadProblem()。此处是ATSP问题，邻接矩阵形式输入;
 * (3)调用 LKH内核 求解：LKHmain().
 *
 * 求解完成后：
 * (1)提取求解结果以返回：函数 OutputTourResult() ;
 * (2)重置 LKH内核 状态：函数 ReSetLKH()。这一步确保当前调用不会对后续调用产生影响.
 */

static char *Copy(char *S);
static void CreateNodes(void);
static int FixEdge(Node * Na, Node * Nb);
static void Read_FIXED_EDGES_SECTION(LKHInput* lkhInput);
static void Read_TOUR_SECTION(LKHInput* lkhInput);
static void Read_EDGE_WEIGHT_SECTION(LKHInput* lkhInput);

void loadDefaultParam(LKHParameters* p){
    p->TimeSpan = 1; 
    p->ScheduleScoreInSecond = 1000;

    p->ProblemFileName = 0;
    p->PiFileName = 0;
    p->InputTourFileName = 0;
    p->OutputTourFileName = 0;
    p->TourFileName = 0;
    p->CandidateFiles = 0;
    p->MergeTourFiles = 0;

    p->AscentCandidates = 50;
    p->BackboneTrials = 0;
    p->Backtracking = 0;
    p->CandidateSetSymmetric = 0;
    p->CandidateSetType = ALPHA;
    p->Crossover = ERXT;
    p->DelaunayPartitioning = 0;
    p->DelaunayPure = 0;
    p->Excess = -1;
    p->ExtraCandidates = 0;
    p->ExtraCandidateSetSymmetric = 0;
    p->ExtraCandidateSetType = QUADRANT;
    p->Gain23Used = 1;
    p->GainCriterionUsed = 1;
    p->GridSize = 1000000.0;
    p->InitialPeriod = -1;
    p->InitialStepSize = 0;
    p->InitialTourAlgorithm = WALK;
    p->InitialTourFraction = 1.0;
    p->KarpPartitioning = 0;
    p->KCenterPartitioning = 0;
    p->KMeansPartitioning = 0;
    p->Kicks = 1;
    p->KickType = 0;
    p->MaxBreadth = INT_MAX;
    p->MaxCandidates = 5;
    p->MaxPopulationSize = 0;
    p->MaxSwaps = -1;
    p->MaxTrials = -1;
    p->MoorePartitioning = 0;
    p->MoveType = 5;
    p->NonsequentialMoveType = -1;
    p->Optimum = MINUS_INFINITY;
    p->PatchingA = 1;
    p->PatchingC = 0;
    p->PatchingAExtended = 0;
    p->PatchingARestricted = 0;
    p->PatchingCExtended = 0;
    p->PatchingCRestricted = 0;
    p->Precision = 100;
    p->POPMUSIC_InitialTour = 0;
    p->POPMUSIC_MaxNeighbors = 5;
    p->POPMUSIC_SampleSize = 10;
    p->POPMUSIC_Solutions = 50;
    p->POPMUSIC_Trials = 1;
    p->Recombination = IPT;
    p->RestrictedSearch = 1;
    p->RohePartitioning = 0;
    p->Runs = 0;
    p->Seed = 1;
    p->SierpinskiPartitioning = 0;
    p->StopAtOptimum = 1;
    p->Subgradient = 1;
    p->SubproblemBorders = 0;
    p->SubproblemsCompressed = 0;
    p->SubproblemSize = 0;
    p->SubsequentMoveType = 0;
    p->SubsequentPatching = 1;
    p->TimeLimit = DBL_MAX;
    p->TotalTimeLimit = DBL_MAX;
    p->TraceLevel = 1;
    p->MaxMatrixDimension = 20006;
}

void ReadParameters(LKHInput* lkhInput)
{
     LKHParameters* lkhParam = lkhInput->lkhParameters;
    // 从结构体获取参数
    TimeSpan = lkhParam->TimeSpan; 
    ScheduleScoreInSecond = lkhParam->ScheduleScoreInSecond;

    ProblemFileName = lkhParam->ProblemFileName;
    PiFileName = lkhParam->PiFileName;
    InputTourFileName = lkhParam->InputTourFileName;
    OutputTourFileName = lkhParam->OutputTourFileName;
    TourFileName = lkhParam->TourFileName;
    CandidateFiles = lkhParam->CandidateFiles;
    MergeTourFiles = lkhParam->MergeTourFiles;
    AscentCandidates = lkhParam->AscentCandidates;
    BackboneTrials = lkhParam->BackboneTrials;
    Backtracking = lkhParam->Backtracking;
    CandidateSetSymmetric = lkhParam->CandidateSetSymmetric;
    CandidateSetType = lkhParam->CandidateSetType;
    Crossover = lkhParam->Crossover;
    DelaunayPartitioning = lkhParam->DelaunayPartitioning;
    DelaunayPure = lkhParam->DelaunayPure;
    Excess = lkhParam->Excess;
    ExtraCandidates = lkhParam->ExtraCandidates;
    ExtraCandidateSetSymmetric = lkhParam->ExtraCandidateSetSymmetric;
    ExtraCandidateSetType = lkhParam->ExtraCandidateSetType;
    Gain23Used = lkhParam->Gain23Used;
    GainCriterionUsed = lkhParam->GainCriterionUsed;
    GridSize = lkhParam->GridSize;
    InitialPeriod = lkhParam->InitialPeriod;
    InitialStepSize = lkhParam->InitialStepSize;
    InitialTourAlgorithm = lkhParam->InitialTourAlgorithm;
    InitialTourFraction = lkhParam->InitialTourFraction;
    KarpPartitioning = lkhParam->KarpPartitioning;
    KCenterPartitioning = lkhParam->KCenterPartitioning;
    KMeansPartitioning = lkhParam->KMeansPartitioning;
    Kicks = lkhParam->Kicks;
    KickType = lkhParam->KickType;
    MaxBreadth = lkhParam->MaxBreadth;
    MaxCandidates = lkhParam->MaxCandidates;
    MaxPopulationSize = lkhParam->MaxPopulationSize;
    MaxSwaps = lkhParam->MaxSwaps;
    MaxTrials = lkhParam->MaxTrials;
    MoorePartitioning = lkhParam->MoorePartitioning;
    MoveType = lkhParam->MoveType;
    NonsequentialMoveType = lkhParam->NonsequentialMoveType;
    Optimum  = lkhParam->Optimum ;
    PatchingA = lkhParam->PatchingA;
    PatchingC = lkhParam->PatchingC;
    PatchingAExtended = lkhParam->PatchingAExtended;
    PatchingARestricted = lkhParam->PatchingARestricted;
    PatchingCExtended = lkhParam->PatchingCExtended;
    PatchingCRestricted = lkhParam->PatchingCRestricted;
    Precision = lkhParam->Precision;
    POPMUSIC_InitialTour = lkhParam->POPMUSIC_InitialTour;
    POPMUSIC_MaxNeighbors = lkhParam->POPMUSIC_MaxNeighbors;
    POPMUSIC_SampleSize = lkhParam->POPMUSIC_SampleSize;
    POPMUSIC_Solutions = lkhParam->POPMUSIC_Solutions;
    POPMUSIC_Trials = lkhParam->POPMUSIC_Trials;
    Recombination = lkhParam->Recombination;
    RestrictedSearch = lkhParam->RestrictedSearch;
    RohePartitioning = lkhParam->RohePartitioning;
    Runs = lkhParam->Runs;
    Seed = lkhParam->Seed;
    SierpinskiPartitioning = lkhParam->SierpinskiPartitioning;
    StopAtOptimum = lkhParam->StopAtOptimum;
    Subgradient = lkhParam->Subgradient;
    SubproblemBorders = lkhParam->SubproblemBorders;
    SubproblemsCompressed = lkhParam->SubproblemsCompressed;
    SubproblemSize = lkhParam->SubproblemSize;
    SubsequentMoveType = lkhParam->SubsequentMoveType;
    SubsequentPatching = lkhParam->SubsequentPatching;
    TimeLimit = lkhParam->TimeLimit;
    TotalTimeLimit = lkhParam->TotalTimeLimit;
    TraceLevel = lkhParam->TraceLevel;
    MaxMatrixDimension = lkhParam->MaxMatrixDimension;

    /* other */
    MergeWithTour =
        Recombination == GPX2 ? MergeWithTourGPX2 :
        Recombination == CLARIST ? MergeWithTourCLARIST :
                                   MergeWithTourIPT;
}

void ReadProblem(LKHInput* lkhInput)
{
    FreeStructures();
    FirstNode = 0;
    C = 0;
    c = 0;

    /* (1) The specification part 定义问题类型、规模和数据输入格式 */
    /* NAME */
    Name = Copy("IOSchedule");
    /* TYPE */
    Type = Copy("ATSP");
    ProblemType = ATSP;
    /* DIMENSION */
    Dimension = lkhInput->matDimension;;
    DimensionSaved = Dimension;
    /* EDGE_WEIGHT_TYPE */
    EdgeWeightType = Copy("EXPLICIT");
    WeightType = EXPLICIT;
    Distance = Distance_EXPLICIT;
    /* EDGE_WEIGHT_FORMAT */
    EdgeWeightFormat = Copy("FULL_MATRIX");
    WeightFormat = FULL_MATRIX;
    /* EDGE_DATA_FORMAT */
    EdgeDataFormat = NULL; 
    /* NODE_COORD_TYPE */
    NodeCoordType = Copy("NO_COORDS");
    CoordType = NO_COORDS;
    /* DISPLAY_DATA_TYPE */
    DisplayDataType = Copy("NO_DISPLAY");
    /* GRID_SIZE */
    GridSize = 1000000.0; // no meaning in here 

    /* (2) The data part 根据上方定义，进行实际数据输入 */
    /* NODE_COORD_SECTION */
    // NULL
    /* EDGE_DATA_SECTION */
    // NULL
    /* FIXED_EDGES_SECTION */
    if(lkhInput->fixEdgeLen != 0)
        Read_FIXED_EDGES_SECTION(lkhInput);
    /* DISPLAY_DATA_SECTION */
    // NULL
    /* TOUR_SECTION */
    if(lkhInput->intialTour != 0)
        Read_TOUR_SECTION(lkhInput);
    /* EDGE_WEIGHT_SECTION */
    Read_EDGE_WEIGHT_SECTION(lkhInput);

    Swaps = 0;

    /* Adjust parameters 调整参数。比如某个参数的默认值依赖于问题类型或输入数据，
       在此前(ReadParameter 函数中)无法确定值是多少，仅是设定为一个默认值标记，
       所以在此处(ReadProblem 函数最后)其默认值才被实际的设定 */
    int i, K;
    if (Seed == 0)
        Seed = (unsigned) (time(0) * (size_t) (&Seed));
    if (Precision == 0)
        Precision = 100;
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (KickType > Dimension / 2)
        KickType = Dimension / 2;
    if (Runs == 0)
        Runs = 10;
    if (MaxCandidates > Dimension - 1)
        MaxCandidates = Dimension - 1;
    if (ExtraCandidates > Dimension - 1)
        ExtraCandidates = Dimension - 1;
    if (SubproblemSize >= Dimension)
        SubproblemSize = Dimension;
    else if (SubproblemSize == 0) {
        if (AscentCandidates > Dimension - 1)
            AscentCandidates = Dimension - 1;
        if (InitialPeriod < 0) {
            InitialPeriod = Dimension / 2;
            if (InitialPeriod < 100)
                InitialPeriod = 100;
        }
        if (Excess < 0)
            Excess = 1.0 / Dimension;
        if (MaxTrials == -1)
            MaxTrials = Dimension;
        HeapMake(Dimension);
    }
    if (POPMUSIC_MaxNeighbors > Dimension - 1)
        POPMUSIC_MaxNeighbors = Dimension - 1;
    if (POPMUSIC_SampleSize > Dimension)
        POPMUSIC_SampleSize = Dimension;
   
    // if (Precision > 1 && (WeightType == EXPLICIT || ProblemType == ATSP)) {
    //     int j, n = ProblemType == ATSP ? Dimension / 2 : Dimension;
    //     for (i = 2; i <= n; i++) {
    //         Node *N = &NodeSet[i];
    //         for (j = 1; j < i; j++)
    //             if (N->C[j] * Precision / Precision != N->C[j])
    //                 eprintf("PRECISION (= %d) is too large", Precision);
    //     }
    // }
    C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
    D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;
    if (SubsequentMoveType == 0)
        SubsequentMoveType = MoveType;
    K = MoveType >= SubsequentMoveType
        || !SubsequentPatching ? MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
        NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1) {
        BestMove = BestSubsequentMove = BestKOptMove;
        if (!SubsequentPatching && SubsequentMoveType <= 5) {
            MoveFunction BestOptMove[] =
                { 0, 0, Best2OptMove, Best3OptMove,
                Best4OptMove, Best5OptMove
            };
            BestSubsequentMove = BestOptMove[SubsequentMoveType];
        }
    } else {
        MoveFunction BestOptMove[] = { 0, 0, Best2OptMove, Best3OptMove,
            Best4OptMove, Best5OptMove
        };
        BestMove = MoveType <= 5 ? BestOptMove[MoveType] : BestKOptMove;
        BestSubsequentMove = SubsequentMoveType <= 5 ?
            BestOptMove[SubsequentMoveType] : BestKOptMove;
    }
    if (TraceLevel >= 1) {
        printff("done\n");
        PrintParameters();
    } 
}

void OutputTourResult(LKHOutput* lkhOutput){
    int j, n = DimensionSaved;
    for (j = 1; j <= n; j++) {
        lkhOutput->tourResult[j - 1] = BestTour[j];
    }
    lkhOutput->tourCost = BestCost;
}

void ReSetLKH(){
    // [TODO 1: 释放 LKH内核 申请(malloc等)的内存, 可以对照源文件 AllocateStructures.c 来做]
    /* tips: 如果在程序结束前没有调用 free，申请的内存将在程序退出时被操作系统回收，
        但这并不算是良好的编程习惯。确保释放内存是管理动态内存的好方法。 */

    // [TODO 2: 重置 LKH内核 的全局变量, 可以对照源文件 LKH.c 来做]
    /* tips: 因为 LKH内核 的运行依赖于全局变量，如果在一个进程内存在多次调用 LKH内核 的话，
        除第一次调用外，后续调用需要重置全局变量，这样才能确保不会对后续调用产生影响 */

    // [TODO 3: 重置 LKH内核 涉及的静态变量, 这些静态变量分布于多个源文件中，需要逐个确认]
    /* tips: 因为 LKH内核 的运行同样依赖于这些静态变量，如果在一个进程内存在多次调用 LKH内核 的话，
        除第一次调用外，后续调用需要重置这些静态变量，这样才能确保不会对后续调用产生影响 */
    /* method: 由于静态变量的作用域仅限于其定义的源文件，无法在其他源文件中直接重新初始化。
        因此，需要使用函数接口：即在定义静态变量的源文件内部，创建一个函数来重新初始化其涉及的静态变量，
        其他源文件就可以通过调用这个函数来实现重新初始化。 */
}

int solveTSP(LKHInput* lkhInput, LKHOutput* lkhOutput){
    /* 参数设定 */
    ReadParameters(lkhInput);
    /* 问题信息设定、数据读入 */
    ReadProblem(lkhInput);
    /* 调用 LKH内核 */
    int ret = LKHmain(lkhInput->scheduleStartTime);
    /* 提取求解结果 */
    OutputTourResult(lkhOutput);
    /* 重置 LKH内核 状态*/
    ReSetLKH();
    return ret;
}

static char *Copy(char *S)
{
    char *Buffer;

    if (!S || strlen(S) == 0)
        return 0;
    Buffer = (char *) malloc(strlen(S) + 1);
    strcpy(Buffer, S);
    return Buffer;
}
static void CreateNodes()
{
    Node *Prev = 0, *N = 0;
    int i;

    if (Dimension <= 0)
        eprintf("DIMENSION is not positive (or not specified)");
    if (ProblemType == ATSP)
        Dimension *= 2;
    NodeSet = (Node *) calloc(Dimension + 1, sizeof(Node));
    for (i = 1; i <= Dimension; i++, Prev = N) {
        N = &NodeSet[i];
        if (i == 1)
            FirstNode = N;
        else
            Link(Prev, N);
        N->Id = i;
        if (MergeTourFiles >= 1)
            N->MergeSuc = (Node **) calloc(MergeTourFiles, sizeof(Node *));
    }
    Link(N, FirstNode);
}
static int FixEdge(Node * Na, Node * Nb)
{
    if (!Na->FixedTo1 || Na->FixedTo1 == Nb)
        Na->FixedTo1 = Nb;
    else if (!Na->FixedTo2 || Na->FixedTo2 == Nb)
        Na->FixedTo2 = Nb;
    else
        return 0;
    if (!Nb->FixedTo1 || Nb->FixedTo1 == Na)
        Nb->FixedTo1 = Na;
    else if (!Nb->FixedTo2 || Nb->FixedTo1 == Na)
        Nb->FixedTo2 = Na;
    else
        return 0;
    return 1;
}
static void Read_FIXED_EDGES_SECTION(LKHInput* lkhInput){
    Node *Ni, *Nj, *N, *NPrev = 0, *NNext;
    int i, j, Count = 0;

    if (!FirstNode)
        CreateNodes();
    int idx = 0;
    while (idx < lkhInput->fixEdgeLen * 2) {
        i = lkhInput->fixEdge[idx++];
        if (i <= 0 || i > DimensionSaved)
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d", i);
        j = lkhInput->fixEdge[idx++];
        if (j <= 0 || j > DimensionSaved)
            eprintf("FIXED_EDGES_SECTION: Node number out of range: %d", j);
        if (i == j)
            eprintf("FIXED_EDGES_SECTION: Illegal edge: %d to %d", i, j);
        Ni = &NodeSet[i];
        Nj = &NodeSet[j + DimensionSaved];
        if (!FixEdge(Ni, Nj))
            eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
        /* Cycle check */
        N = Ni;
        Count = 0;
        do {
            NNext = N->FixedTo1 != NPrev ? N->FixedTo1 : N->FixedTo2;
            NPrev = N;
            Count++;
        } while ((N = NNext) && N != Ni);
        if (N == Ni && Count != Dimension)
            eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
    }
}
static void Read_TOUR_SECTION(LKHInput* lkhInput) {
    Node* First = 0, * Last = 0, * N, * Na;
    int i, k;
    
    if (TraceLevel >= 1) {
        printff("Reading TOUR...");
    }
    if (!FirstNode)
        CreateNodes();
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);

    int idx = 0;
    i = lkhInput->intialTour[idx];
    for (k = 0; k <= DimensionSaved ; k++) {
        if (i <= 0 || i > DimensionSaved)
            eprintf("(TOUR_SECTION) Node number out of range: %d",i);
        N = &NodeSet[i];
        if (N->V == 1 && k != DimensionSaved)
            eprintf("(TOUR_SECTION) Node number occurs twice: %d",N->Id);
        N->V = 1;
        if (k == 0)
            First = Last = N;
        else {
            Na = N + DimensionSaved;
            Na->V = 1;
            /* 划分子问题需要 */
            if (!Na)
                (Last->SubproblemSuc = N)->SubproblemPred = Last;
            else {
                (Last->SubproblemSuc = Na)->SubproblemPred = Last;
                (Na->SubproblemSuc = N)->SubproblemPred = Na;
            }
            /* 初始解需要 */
            if (!Na)
                Last->InitialSuc = N;
            else {
                Last->InitialSuc = Na;
                Na->InitialSuc = N;
            }

            Last = N;
        }
        if (k < DimensionSaved) {
            i = lkhInput->intialTour[++idx];
        }
        if (k == DimensionSaved - 1)
            i = First->Id;
    }
  
    N = FirstNode;
    do {    
        if (!N->V)
            printf("TOUR_SECTION: Node is missing: %d", N->Id);
    } while ((N = N->Suc) != FirstNode);
        do {
            if (N->FixedTo1 &&
                N->SubproblemPred != N->FixedTo1
                && N->SubproblemSuc != N->FixedTo1)
                eprintf("Fixed edge (%d, %d) "
                    "does not belong to subproblem tour", N->Id,
                    N->FixedTo1->Id);
            if (N->FixedTo2 && N->SubproblemPred != N->FixedTo2
                && N->SubproblemSuc != N->FixedTo2)
                eprintf("Fixed edge (%d, %d) "
                    "does not belong to subproblem tour", N->Id,
                    N->FixedTo2->Id);
        } while ((N = N->Suc) != FirstNode);
    if (TraceLevel >= 1)
        printff("done\n");
}
static void Read_EDGE_WEIGHT_SECTION(LKHInput* lkhInput)
{
    Node *Ni, *Nj;
    int i, j, n, W;

    if (!FirstNode)
        CreateNodes(); // 创建结点数组，是一个头尾相连的双链表数据结构，FirstNode指向第一个结点

    /* ProblemType == ATSP && WeightFormat == FULL_MATRIX */
    // 为矩阵申请空间
    n = DimensionSaved;
    CostMatrix = (int *) calloc((size_t) n * n, sizeof(int));
    for (Ni = FirstNode; Ni->Id <= n; Ni = Ni->Suc)
        Ni->C = &CostMatrix[(size_t) (Ni->Id - 1) * n] - 1;
    // 写入矩阵中
    for (i = 1; i <= n; i++) {
        Ni = &NodeSet[i];
        for (j = 1; j <= n; j++) {
            W = lkhInput->adjMat[i - 1][j - 1]; // 对应位置
            if (j != i && W > INT_MAX / 2 / Precision)
                eprintf("EDGE_WEIGHT_SECTION: "
                        "Weight %d > INT_MAX / 2 / PRECISION", W);
            Ni->C[j] = W;
            if (i != j && W > M)
                M = W;
        }
        Nj = &NodeSet[i + n];
        FixEdge(Ni, Nj);
    }
    Distance = Distance_ATSP;
    WeightType = -1;
}
//################################### LKHInterface.c end ###################################


//################################### Activate.c begin ###################################
#include "LKH.h"

/*
 * If a node is made "active", attempts are made to find an improving 
 * move with the node as anchor node, t1.
 *
 * The Active function makes a node active by inserting it into a 
 * queue of active nodes. FirstActive denotes the front node of 
 * the queue. LastActive denotes the rear. 
 *
 * The queue is implemented as a circular list in which the Next field 
 * of each Node references the successor node. 
 *
 * A node is member of the queue iff its Next != 0. The function has no 
 * effect if the node is already in the queue. 
 *
 * The function is called from the StoreTour function.  
 */

void Activate(Node * N)
{
    if (N->Next != 0)
        return;
    if (FixedOrCommon(N, N->Pred) && FixedOrCommon(N, N->Suc))
        return;
    if (FirstActive == 0)
        FirstActive = LastActive = N;
    else
        LastActive = LastActive->Next = N;
    LastActive->Next = FirstActive;
}
//################################### Activate.c end ###################################

//################################### AddCandidate.c begin ###################################
#include "LKH.h"

/*
 * The AddCandidate function adds a given edge (From, To) to the set
 * of candidate edges associated with the node From. The cost and 
 * alpha-value of the edge are passed as parameters to the function.
 *
 * The function has no effect if the edge is already in the candidate
 * set.
 *
 * If the edge was added, the function returns 1; otherwise 0.
 *    
 * The function is called from the functions CreateDelaunaySet and
 * OrderCandidateSet.
 */

int AddCandidate(Node * From, Node * To, int Cost, int Alpha)
{
    int Count;
    Candidate *NFrom;

    if (From->Subproblem != FirstNode->Subproblem ||
        To->Subproblem != FirstNode->Subproblem ||
        Cost == INT_MAX)
        return 0;
    if (From->CandidateSet == 0)
        From->CandidateSet = (Candidate *) calloc(3, sizeof(Candidate));
    if (From == To || To->Subproblem != FirstNode->Subproblem ||
        !IsPossibleCandidate(From, To))
        return 0;
    Count = 0;
    for (NFrom = From->CandidateSet; NFrom->To && NFrom->To != To; NFrom++)
        Count++;
    if (NFrom->To) {
        if (NFrom->Alpha == INT_MAX)
            NFrom->Alpha = Alpha;
        return 0;
    }
    NFrom->Cost = Cost;
    NFrom->Alpha = Alpha;
    NFrom->To = To;
    From->CandidateSet =
        (Candidate *) realloc(From->CandidateSet,
                              (Count + 2) * sizeof(Candidate));
    From->CandidateSet[Count + 1].To = 0;
    return 1;
}
//################################### AddCandidate.c end ###################################

//################################### AddTourCandidates.c begin ###################################
#include "LKH.h"

/*
 * The AddTourCandidates function extends the candidate set with tour
 * edges given in the tour files.
 *   
 * The function is called from GenerateCandidateSet and OrderCandidateSet.  
*/

void AddTourCandidates()
{
    Node *Na, *Nb;
    int i, d, Subproblem = FirstNode->Subproblem;

    /* Add fixed edges */
    Na = FirstNode;
    do {
        if (Na->FixedTo1)
            AddCandidate(Na, Na->FixedTo1, D(Na, Na->FixedTo1), 0);
        if (Na->FixedTo2)
            AddCandidate(Na, Na->FixedTo2, D(Na, Na->FixedTo2), 0);
    }
    while ((Na = Na->Suc) != FirstNode);

    /* Add MERGE_TOUR_FILE edges */
    for (i = 0; i < MergeTourFiles; i++) {
        Na = FirstNode;
        do {
            Nb = Na->MergeSuc[i];
            if (!Nb)
                break;
            if (Na->Subproblem == Subproblem &&
                Nb->Subproblem == Subproblem) {
                d = D(Na, Nb);
                AddCandidate(Na, Nb, d, 1);
                AddCandidate(Nb, Na, d, 1);
            }
        }
        while ((Na = Nb) != FirstNode);
    }

    /* Add INITIAL_TOUR_FILE edges */
    Na = FirstNode;
    do {
        Nb = Na->InitialSuc;
        if (!Nb)
            break;
        if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
            d = D(Na, Nb);
            AddCandidate(Na, Nb, d, 1);
            AddCandidate(Nb, Na, d, 1);
        }
    }
    while ((Na = Nb) != FirstNode);

    /* Add INPUT_TOUR_FILE edges */
    Na = FirstNode;
    do {
        Nb = Na->InputSuc;
        if (!Nb)
            break;
        if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
            d = D(Na, Nb);
            AddCandidate(Na, Nb, d, 1);
            AddCandidate(Nb, Na, d, 1);
        }
    }
    while ((Na = Nb) != FirstNode);

    /* Add SUBPROBLEM_TOUR_FILE edges */
    Na = FirstNode;
    do {
        Nb = Na->SubproblemSuc;
        if (!Nb)
            break;
        if (Na->Subproblem == Subproblem && Nb->Subproblem == Subproblem) {
            d = D(Na, Nb);
            AddCandidate(Na, Nb, d, 1);
            AddCandidate(Nb, Na, d, 1);
        }
    } while ((Na = Nb) != FirstNode);
}
//################################### AddTourCandidates.c end ###################################

//################################### AdjustCandidateSet.c begin ###################################
#include "LKH.h"

/*
 * When a trial has produced a better tour, the set of candidate edges is
 * adjusted as follows:
 *
 *   (1) The set is extended with tour edges not present in the
 *       current set;
 *   (2) Precedence is given to those edges that are common to the two
 *       currently best tours.
 *
 * The AdjustCandidateSet function adjusts for each node its table of 
 * candidate edges. A new candidate edge is added by extending the table 
 * and inserting the edge as its last ordinary element (disregarding the 
 * dummy edge). The Alpha field of the new candidate edge is set to 
 * INT_MAX. Edges that belong to the best tour as well as the next best 
 * tour are moved to the start of the table.                         
 */

void AdjustCandidateSet()
{
    Candidate *NFrom, *NN, Temp;
    Node *From = FirstNode, *To;

    /* Extend and reorder candidate sets */
    do {
        if (!From->CandidateSet)
            From->CandidateSet = (Candidate *) calloc(3, sizeof(Candidate));
        /* Extend */
        for (To = From->Pred; To; To = To == From->Pred ? From->Suc : 0) {
            int Count = 0;
            if ((ProblemType == HCP || ProblemType == HPP) &&
                !IsBackboneCandidate(From, To))
                continue;
            for (NFrom = From->CandidateSet; NFrom->To && NFrom->To != To;
                 NFrom++)
                Count++;
            if (!NFrom->To) {
                /* Add new candidate edge */
                NFrom->Cost = C(From, To);
                NFrom->To = To;
                NFrom->Alpha = INT_MAX;
                From->CandidateSet =
                   (Candidate *) realloc(From->CandidateSet,
                                         (Count + 2) * sizeof(Candidate));
                From->CandidateSet[Count + 1].To = 0;
            }
        }
        /* Reorder */
        for (NFrom = From->CandidateSet + 1; (To = NFrom->To); NFrom++)
            if (InBestTour(From, To) && InNextBestTour(From, To)) {
                /* Move the edge to the start of the candidate table */
                Temp = *NFrom;
                for (NN = NFrom - 1; NN >= From->CandidateSet; NN--)
                    *(NN + 1) = *NN;
                *(NN + 1) = Temp;
            }
    }
    while ((From = From->Suc) != FirstNode);
}
//################################### AdjustCandidateSet.c end ###################################

//################################### AllocateStructures.c begin ###################################
#include "Segment.h"
#include "LKH.h"
#include "Heap.h"
#include "Sequence.h"

/*      
 * The AllocateStructures function allocates all necessary 
 * structures except nodes and candidates.
 */

#define Free(s) { free(s); s = 0; }

void AllocateStructures()
{
    int i, K;

    Free(Heap);
    Free(BestTour);
    Free(BetterTour);
    Free(HTable);
    Free(Rand);
    Free(CacheSig);
    Free(CacheVal);
    Free(T);
    Free(G);
    Free(t);
    Free(p);
    Free(q);
    Free(SwapStack);
    Free(tSaved);

    HeapMake(Dimension);
    BestTour = (int *) calloc(1 + Dimension, sizeof(int));
    BetterTour = (int *) calloc(1 + Dimension, sizeof(int));
    HTable = (HashTable *) malloc(sizeof(HashTable));
    HashInitialize((HashTable *) HTable);
    SRandom(Seed);
    Rand = (unsigned *) malloc((Dimension + 1) * sizeof(unsigned));
    for (i = 1; i <= Dimension; i++)
        Rand[i] = Random();
    SRandom(Seed);
    if (WeightType != EXPLICIT) {
        for (i = 0; (1 << i) < (Dimension << 1); i++);
        i = 1 << i;
        CacheSig = (int *) calloc(i, sizeof(int));
        CacheVal = (int *) calloc(i, sizeof(int));
        CacheMask = i - 1;
    }
    AllocateSegments();
    K = MoveType;
    if (SubsequentMoveType > K)
        K = SubsequentMoveType;
    T = (Node **) malloc((1 + 2 * K) * sizeof(Node *));
    G = (GainType *) malloc(2 * K * sizeof(GainType));
    t = (Node **) malloc(6 * K * sizeof(Node *));
    tSaved = (Node **) malloc((1 + 2 * K) * sizeof(Node *));
    p = (int *) malloc(6 * K * sizeof(int));
    q = (int *) malloc(6 * K * sizeof(int));
    incl = (int *) malloc(6 * K * sizeof(int));
    cycle = (int *) malloc(6 * K * sizeof(int));
    SwapStack =
        (SwapRecord *) malloc((MaxSwaps + 6 * K) * sizeof(SwapRecord));
}

/*      
 * The AllocateSegments function allocates the segments of the two-level tree.
 */

void AllocateSegments()
{
    Segment *S = 0, *SPrev;
    SSegment *SS = 0, *SSPrev;
    int i;

    FreeSegments();
#ifdef THREE_LEVEL_TREE
    GroupSize = (int) pow((double) Dimension, 1.0 / 3.0);
#elif defined TWO_LEVEL_TREE
    GroupSize = (int) sqrt((double) Dimension);
#else
    GroupSize = Dimension;
#endif
    Groups = 0;
    for (i = Dimension, SPrev = 0; i > 0; i -= GroupSize, SPrev = S) {
        S = (Segment *) malloc(sizeof(Segment));
        S->Rank = ++Groups;
        if (!SPrev)
            FirstSegment = S;
        else
            SLink(SPrev, S);
    }
    SLink(S, FirstSegment);
#ifdef THREE_LEVEL_TREE
    SGroupSize = sqrt((double) Groups);
#else
    SGroupSize = Dimension;
#endif
    SGroups = 0;
    for (i = Groups, SSPrev = 0; i > 0; i -= SGroupSize, SSPrev = SS) {
        SS = (SSegment *) malloc(sizeof(SSegment));
        SS->Rank = ++SGroups;
        if (!SSPrev)
            FirstSSegment = SS;
        else
            SLink(SSPrev, SS);
    }
    SLink(SS, FirstSSegment);
}
//################################### AllocateStructures.c end ###################################

//################################### Ascent.c begin ###################################
#include "LKH.h"

/* 
 * The Ascent function computes a lower bound on the optimal tour length 
 * using subgradient optimization. The function also transforms the original 
 * problem into a problem in which the Alpha-values reflect the likelihood 
 * of edges being optimal.
 *
 * The function attempts to find penalties (Pi-values) that maximizes the 
 * lower bound L(T(Pi)) - 2*PiSum, where L(T(Pi)) denotes the length of the 
 * minimum spanning 1-tree computed from the transformed distances, and PiSum 
 * denotes the sum of Pi-values. If C(i,j) denotes the length of an edge 
 * (i,j), then the transformed distance D(i,j) of an edge is 
 * C(i,j) + Pi(i) + Pi(j).
 *
 * The Minimum1TreeCost function is used to compute the cost of a minimum 
 * 1-tree.The Generatecandidates function is called in order to generate 
 * candidate sets. Minimum 1-trees are then computed in the corresponding 
 * sparse graph.         
 */

GainType Ascent()
{
    Node *t;
    GainType BestW, W, W0, Alpha, MaxAlpha;
    int T, Period, P, InitialPhase, BestNorm;

  Start:
    /* Initialize Pi and BestPi */
    t = FirstNode;
    do
        t->Pi = t->BestPi = 0;
    while ((t = t->Suc) != FirstNode);
    
    if (CandidateSetType == POPMUSIC && !FirstNode->CandidateSet)
        Create_POPMUSIC_CandidateSet(AscentCandidates);
    else if (MaxCandidates == 0) {
        AddTourCandidates();
    }

    /* Compute the cost of a minimum 1-tree */
    W = Minimum1TreeCost(CandidateSetType == POPMUSIC ||
                         MaxCandidates == 0);

    /* Return this cost 
       if either
       (1) subgradient optimization is not wanted, or
       (2) the norm of the tree (its deviation from a tour) is zero
       (in that case the true optimum has been found).
     */
    if (!Subgradient || !Norm)
        return W;

    if (MaxCandidates > 0) {
        /* Generate symmetric candididate sets for all nodes */
        MaxAlpha = INT_MAX;
        if (Optimum != MINUS_INFINITY
            && (Alpha = Optimum * Precision - W) >= 0)
            MaxAlpha = Alpha;
        if (CandidateSetType != POPMUSIC)
            GenerateCandidates(AscentCandidates, MaxAlpha, 1);
        else {
            OrderCandidateSet(AscentCandidates, MaxAlpha, 1);
            W = Minimum1TreeCost(1);
            if (!Norm || W / Precision == Optimum)
                return W;
        }
    }
    if (TraceLevel >= 2) {
        CandidateReport();
        printff("Subgradient optimization ...\n");
    }

    /* Set LastV of every node to V (the node's degree in the 1-tree) */
    t = FirstNode;
    do
        t->LastV = t->V;
    while ((t = t->Suc) != FirstNode);

    BestW = W0 = W;
    BestNorm = Norm;
    InitialPhase = 1;
    /* Perform subradient optimization with decreasing period length 
       and decreasing step size */
    for (Period = InitialPeriod, T = InitialStepSize * Precision;
         Period > 0 && T > 0 && Norm != 0; Period /= 2, T /= 2) {
        /* Period and step size are halved at each iteration */
        if (TraceLevel >= 2)
            printff
                ("  T = %d, Period = %d, BestW = %0.1f, BestNorm = %d\n",
                 T, Period, (double) BestW / Precision, BestNorm);
        for (P = 1; T && P <= Period && Norm != 0; P++) {
            /* Adjust the Pi-values */
            t = FirstNode;
            do {
                if (t->V != 0) {
                    t->Pi += T * (7 * t->V + 3 * t->LastV) / 10;
                    if (t->Pi > INT_MAX / 10)
                        t->Pi = INT_MAX / 10;
                    else if (t->Pi < INT_MIN / 10)
                        t->Pi = INT_MIN / 10;
                }
                t->LastV = t->V;
            }
            while ((t = t->Suc) != FirstNode);
            /* Compute a minimum 1-tree in the sparse graph */
            W = Minimum1TreeCost(1);
            /* Test if an improvement has been found */
            if (W > BestW || (W == BestW && Norm < BestNorm)) {
                /* If the lower bound becomes greater than twice its
                   initial value it is taken as a sign that the graph might be
                   too sparse */
                if (W - W0 > (W0 >= 0 ? W0 : -W0) && AscentCandidates > 0
                    && AscentCandidates < Dimension) {
                    W = Minimum1TreeCost(CandidateSetType == DELAUNAY ||
                                         CandidateSetType == POPMUSIC ||
                                         MaxCandidates == 0);
                    if (W < W0) {
                        /* Double the number of candidate edges 
                           and start all over again */
                        if (TraceLevel >= 2)
                            printff("Warning: AscentCandidates doubled\n");
                        if ((AscentCandidates *= 2) > Dimension)
                            AscentCandidates = Dimension;
                        goto Start;
                    }
                    W0 = W;
                }
                BestW = W;
                BestNorm = Norm;
                /* Update the BestPi-values */
                t = FirstNode;
                do
                    t->BestPi = t->Pi;
                while ((t = t->Suc) != FirstNode);
                if (TraceLevel >= 2)
                    printff
                        ("* T = %d, Period = %d, P = %d, "
                         "BestW = %0.1f, BestNorm = %d\n",
                         T, Period, P, (double) BestW / Precision,
                         BestNorm);
                /* If in the initial phase, the step size is doubled */
                if (InitialPhase && T * sqrt((double) Norm) > 0)
                    T *= 2;
                /* If the improvement was found at the last iteration of the 
                   current period, then double the period */
                if (CandidateSetType != POPMUSIC &&
                    P == Period && (Period *= 2) > InitialPeriod)
                    Period = InitialPeriod;
            } else {
                if (TraceLevel >= 3)
                    printff
                        ("  T = %d, Period = %d, P = %d, W = %0.1f, Norm = %d\n",
                         T, Period, P, (double) W / Precision, Norm);
                if (InitialPhase && P > Period / 2) {
                    /* Conclude the initial phase */
                    InitialPhase = 0;
                    P = 0;
                    T = 3 * T / 4;
                }
            }
        }
    }

    t = FirstNode;
    do {
        t->Pi = t->BestPi;
        t->BestPi = 0;
    } while ((t = t->Suc) != FirstNode);

    /* Compute a minimum 1-tree */
    W = BestW = Minimum1TreeCost(CandidateSetType == POPMUSIC ||
                                 MaxCandidates == 0);

    if (MaxCandidates > 0 && CandidateSetType != POPMUSIC) {
        FreeCandidateSets();
    } else {
        Candidate *Nt;
        t = FirstNode;
        do {
            for (Nt = t->CandidateSet; Nt && Nt->To; Nt++)
                Nt->Cost += t->Pi + Nt->To->Pi;
        }
        while ((t = t->Suc) != FirstNode);
    }
    if (TraceLevel >= 2)
        printff("Ascent: BestW = %0.1f, Norm = %d\n",
                (double) BestW / Precision, Norm);
    return W;
}
//################################### Ascent.c end ###################################

//################################### Best2OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Best2OptMove function makes sequential edge exchanges. If possible, it 
 * makes a 2-opt move that improves the tour. Otherwise, it makes the most 
 * promising 2-opt move that fulfils the positive gain criterion. To prevent 
 * an infinity chain of moves the last edge in a 2-opt move must not previously 
 * have been included in the chain. 
 *
 * The edge (t1,t2) is the first edge to be exchanged. G0 is a pointer to the 
 * accumulated gain.
 *
 * In case a 2-opt move is found that improves the tour, the improvement of 
 * the cost is made available to the caller through the parameter Gain. 
 * If *Gain > 0, an improvement of the current tour has been found. In this
 * case the function returns 0.
 *
 * Otherwise, the best 2-opt move is made, and a pointer to the node that was 
 * connected to t1 (in order to close the tour) is returned. The new 
 * accumulated gain is made available to the caller through the parameter G0. 
 *
 * The function is called from the LinKernighan function. 
 */

Node *Best2OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
{
    Node *t3, *t4, *T3 = 0, *T4 = 0;
    Candidate *Nt2;
    GainType G1, G2, BestG2 = MINUS_INFINITY;
    int Breadth2 = 0;

    if (ProblemType == ATSP)
        return 0;
    if (SUC(t1) != t2)
        Reversed ^= 1;

    /* 
     * Determine (T3,T4) = (t3,t4)
     * such that 
     *
     *     G4 = *G0 - C(t2,T3) + C(T3,T4)
     *
     * is maximum (= BestG2), and (T3,T4) has not previously been included.
     * If during this process a legal move with *Gain > 0 is found, then make
     * the move and exit Best2OptMove immediately 
     */

    /* Choose (t2,t3) as a candidate edge emanating from t2 */
    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        if (t3 == t2->Pred || t3 == t2->Suc ||
            ((G1 = *G0 - Nt2->Cost) <= 0 && GainCriterionUsed &&
             ProblemType != HCP && ProblemType != HPP))
            continue;
        /* Choose t4 (only one choice gives a closed tour) */
        t4 = PRED(t3);
        if (FixedOrCommon(t3, t4))
            continue;
        G2 = G1 + C(t3, t4);
        if (!Forbidden(t4, t1) &&
            (!c || G2 - c(t4, t1) > 0) && (*Gain = G2 - C(t4, t1)) > 0) {
            Swap1(t1, t2, t3);
            return 0;
        }
        if (++Breadth2 > MaxBreadth)
            break;
        if (GainCriterionUsed && G2 - Precision < t4->Cost)
            continue;
        if (!Backtracking || Swaps > 0) {
            if ((G2 > BestG2 ||
                 (G2 == BestG2 && !Near(t3, t4) &&
                  Near(T3, T4))) &&
                Swaps < MaxSwaps &&
                Excludable(t3, t4) && !InInputTour(t3, t4)) {
                T3 = t3;
                T4 = t4;
                BestG2 = G2;
            }
        } else if (MaxSwaps > 0) {
            GainType G = G2;
            Node *t = t4;
            Make2OptMove(t1, t2, t3, t4);
            Exclude(t1, t2);
            Exclude(t3, t4);
            while ((t = BestSubsequentMove(t1, t, &G, Gain)));
            if (*Gain > 0)
                return 0;
            RestoreTour();
            if (t2 != SUC(t1))
                Reversed ^= 1;
        }
    }
    *Gain = 0;
    if (T4) {
        /* Make the best 2-opt move */
        Swap1(t1, t2, T3);
        Exclude(t1, t2);
        Exclude(T3, T4);
        *G0 = BestG2;
    }
    return T4;
}

/*
   Below is shown the use of the variable X4 to discriminate between 
   the 2 cases considered by the algorithm. 

   The notation

        ab-

   is used for a subtour that starts with the edge (ta,tb). For example 
   the tour 

        12-43-

   contains the edges (t1,t2) and (t4,t3), in that order. 

   X4 = 1:
       12-43-
   X4 = 2:
       12-34-
*/
//################################### Best2OptMove.c end ###################################

//################################### Best3OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Best3OptMove function makes sequential edge exchanges. If possible, it 
 * makes an  r-opt move (r <= 3) that improves the tour. Otherwise, it makes 
 * the most promising 3-opt move that fulfils the positive gain criterion. 
 * To prevent an infinity chain of moves the last edge in a 3-opt move must 
 * not previously have been included in the chain. 
 *
 * The edge (t1,t2) is the first edge to be exchanged. G0 is a pointer to the 
 * accumulated gain.
 *
 * In case a r-opt move is found that improves the tour, the improvement of 
 * the cost is made available to the caller through the parameter Gain. 
 * If *Gain > 0, an improvement of the current tour has been found. In this
 * case the function returns 0.
 *
 * Otherwise, the best 3-opt move is made, and a pointer to the node that was 
 * connected to t1 (in order to close the tour) is returned. The new 
 * accumulated gain is made available to the caller through the parameter G0. 
 *
 * The function is called from the LinKernighan function. 
 */

/* 
   The algorithm splits the set of possible moves up into a number disjoint 
   subsets (called "cases"). When t1, t2, ..., t6 has been chosen, Case6 is 
   used to discriminate between 8 cases.
 
   A description of the cases is given after the code.   
*/

Node *Best3OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
{
    Node *t3, *t4, *t5, *t6, *T3 = 0, *T4 = 0, *T5 = 0, *T6 = 0;
    Candidate *Nt2, *Nt4;
    GainType G1, G2, G3, G4, BestG4 = MINUS_INFINITY;
    int Case6, BestCase6 = 0, X4, X6;
    int Breadth2 = 0, Breadth4;

    if (SUC(t1) != t2)
        Reversed ^= 1;

    /* 
     * Determine (T3,T4,T5,T6) = (t3,t4,t5,t6)
     * such that 
     *
     *     G4 = *G0 - C(t2,T3) + C(T3,T4) 
     *              - C(T4,T5) + C(T5,T6)
     *
     * is maximum (= BestG4), and (T5,T6) has not previously been included.
     * If during this process a legal move with *Gain > 0 is found, then make
     * the move and exit Best3OptMove immediately. 
     */

    /* Choose (t2,t3) as a candidate edge emanating from t2 */
    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        if (t3 == t2->Pred || t3 == t2->Suc ||
            ((G1 = *G0 - Nt2->Cost) <= 0 && GainCriterionUsed &&
             ProblemType != HCP && ProblemType != HPP))
            continue;
        if (++Breadth2 > MaxBreadth)
            break;
        /* Choose t4 as one of t3's two neighbors on the tour */
        for (X4 = ProblemType == ATSP ? 2 : 1; X4 <= 2; X4++) {
            t4 = X4 == 1 ? PRED(t3) : SUC(t3);
            if (FixedOrCommon(t3, t4))
                continue;
            G2 = G1 + C(t3, t4);
            if (X4 == 1 && !Forbidden(t4, t1) &&
                (!c || G2 - c(t4, t1) > 0) && (*Gain = G2 - C(t4, t1)) > 0)
            {
                Swap1(t1, t2, t3);
                return 0;
            }
            if (Backtracking && !Excludable(t3, t4))
                continue;
            Breadth4 = 0;
            /* Choose (t4,t5) as a candidate edge emanating from t4 */
            for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {
                if (t5 == t4->Pred || t5 == t4->Suc ||
                    ((G3 = G2 - Nt4->Cost) <= 0 && GainCriterionUsed &&
                     ProblemType != HCP && ProblemType != HPP) ||
                    (X4 == 2 && !BETWEEN(t2, t5, t3)))
                    continue;
                if (++Breadth4 > MaxBreadth)
                    break;
                /* Choose t6 as one of t5's two neighbors on the tour */
                for (X6 = 1; X6 <= X4; X6++) {
                    if (X4 == 1) {
                        Case6 = 1 + !BETWEEN(t2, t5, t4);
                        t6 = Case6 == 1 ? SUC(t5) : PRED(t5);
                    } else {
                        Case6 = 4 + X6;
                        t6 = X6 == 1 ? SUC(t5) : PRED(t5);
                        if (t6 == t1)
                            continue;
                    }
                    if (FixedOrCommon(t5, t6))
                        continue;
                    G4 = G3 + C(t5, t6);
                    if (!Forbidden(t6, t1) &&
                        (!c || G4 - c(t6, t1) > 0) &&
                        (*Gain = G4 - C(t6, t1)) > 0) {
                        Make3OptMove(t1, t2, t3, t4, t5, t6, Case6);
                        return 0;
                    }
                    if (GainCriterionUsed && G4 - Precision < t6->Cost)
                        continue;
                    if (!Backtracking || Swaps > 0) {
                        if ((G4 > BestG4 ||
                             (G4 == BestG4 && !Near(t5, t6) &&
                              Near(T5, T6))) &&
                            Swaps < MaxSwaps &&
                            Excludable(t5, t6) && !InInputTour(t5, t6)) {
                            /* Ignore the move if the gain does not vary */
                            if (RestrictedSearch &&
                                ProblemType != HCP &&
                                ProblemType != HPP &&
                                G2 - t4->Pi == G4 - t6->Pi &&
                                G3 + t5->Pi == G1 + t3->Pi)
                                continue;
                            T3 = t3;
                            T4 = t4;
                            T5 = t5;
                            T6 = t6;
                            BestCase6 = Case6;
                            BestG4 = G4;
                        }
                    } else if (MaxSwaps > 0) {
                        GainType G = G4;
                        Node *t = t6;
                        Make3OptMove(t1, t2, t3, t4, t5, t6, Case6);
                        Exclude(t1, t2);
                        Exclude(t3, t4);
                        Exclude(t5, t6);
                        while ((t = BestSubsequentMove(t1, t, &G, Gain)));
                        if (*Gain > 0)
                            return 0;
                        RestoreTour();
                        if (t2 != SUC(t1))
                            Reversed ^= 1;
                    }
                }
            }
        }
    }
    *Gain = 0;
    if (T6) {
        /* Make the best 3-opt move */
        Make3OptMove(t1, t2, T3, T4, T5, T6, BestCase6);
        Exclude(t1, t2);
        Exclude(T3, T4);
        Exclude(T5, T6);
        *G0 = BestG4;
    }
    return T6;
}

/*
   Below is shown the use of the variables X4 and Case6 to discriminate between 
   the 4 cases considered by the algorithm. 

   The notation

        ab-

   is used for a subtour that starts with the edge (ta,tb). For example 
   the tour 

        12-43-

   contains the edges (t1,t2) and (t4,t3), in that order. 

   X4 = 1:
       12-43-
       Case6 = 1: 
           12-56-43-
       Case6 = 2:   
           12-43-65-
   X4 = 2:
       12-34-
       Case6 = 5: 
           12-56-34-
       Case6 = 6: 
           12-65-34-
*/
//################################### Best3OptMove.c end ###################################

//################################### Best4OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Best4OptMove function makes sequential edge exchanges. If possible, it 
 * makes an r-opt move (r <= 4) that improves the tour. Otherwise, it makes 
 * the most promising 4-opt move that fulfils the positive gain criterion. 
 * To prevent an infinity chain of moves the last edge in a 4-opt move must 
 * not previously have been included in the chain. 
 *
 * The edge (t1,t2) is the first edge to be exchanged. G0 is a pointer to the 
 * accumulated gain.
 *
 * In case a r-opt move is found that improves the tour, the improvement of 
 * the cost is made available to the caller through the parameter Gain. 
 * If *Gain > 0, an improvement of the current tour has been found. In this
 * case the function returns 0.
 *
 * Otherwise, the best 4-opt move is made, and a pointer to the node that was 
 * connected to t1 (in order to close the tour) is returned. The new 
 * accumulated gain is made available to the caller through the parameter G0. 
 *
 * The function is called from the LinKernighan function. 
 */

/* 
   The algorithm splits the set of possible moves up into a number disjoint 
   subsets (called "cases"). When t1, t2, ..., t6 has been chosen, Case6 is 
   used to discriminate between 8 cases. When t1, t2, ..., t8 has been chosen, 
   Case8 is used to discriminate between 16 cases. 

   A description of the cases is given after the code.   
*/

Node *Best4OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
{
    Node *t3, *t4, *t5, *t6 = 0, *t7, *t8 = 0,
        *T3 = 0, *T4 = 0, *T5 = 0, *T6 = 0, *T7 = 0, *T8 = 0;
    Candidate *Nt2, *Nt4, *Nt6;
    GainType G1, G2, G3, G4, G5, G6, BestG6 = MINUS_INFINITY;
    int Case6 = 0, Case8 = 0, BestCase8 = 0, X4, X6, X8;
    int Breadth2 = 0, Breadth4, Breadth6;

    if (SUC(t1) != t2)
        Reversed ^= 1;

    /* 
     * Determine (T3,T4,T5,T6,T7,T8) = (t3,t4,t5,t6,t7,t8)
     * such that
     *
     *     G8 = *G0 - C(t2,T3) + C(T3,T4) 
     *              - C(T4,T5) + C(T5,T6)
     *              - C(T6,T7) + C(T7,T8)
     *
     * is maximum (= BestG6), and (T7,T8) has not previously been included.
     * If during this process a legal move with *Gain > 0 is found, then make
     * the move and exit Best4OptMove immediately. 
     */

    /* Choose (t2,t3) as a candidate edge emanating from t2 */
    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        if (t3 == t2->Pred || t3 == t2->Suc ||
            ((G1 = *G0 - Nt2->Cost) <= 0 && GainCriterionUsed &&
             ProblemType != HCP && ProblemType != HPP))
            continue;
        if (++Breadth2 > MaxBreadth)
            break;
        /* Choose t4 as one of t3's two neighbors on the tour */
        for (X4 = ProblemType == ATSP ? 2 : 1; X4 <= 2; X4++) {
            t4 = X4 == 1 ? PRED(t3) : SUC(t3);
            if (FixedOrCommon(t3, t4))
                continue;
            G2 = G1 + C(t3, t4);
            if (X4 == 1 && !Forbidden(t4, t1) &&
                (!c || G2 - c(t4, t1) > 0) && (*Gain = G2 - C(t4, t1)) > 0)
            {
                Swap1(t1, t2, t3);
                return 0;
            }
            if (Backtracking && !Excludable(t3, t4))
                continue;
            Breadth4 = 0;
            /* Choose (t4,t5) as a candidate edge emanating from t4 */
            for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {
                if (t5 == t4->Pred || t5 == t4->Suc ||
                    ((G3 = G2 - Nt4->Cost) <= 0 && GainCriterionUsed &&
                     ProblemType != HCP && ProblemType != HPP))
                    continue;
                if (++Breadth4 > MaxBreadth)
                    break;
                /* Choose t6 as one of t5's two neighbors on the tour */
                for (X6 = 1; X6 <= 2; X6++) {
                    if (X4 == 1) {
                        if (X6 == 1) {
                            Case6 = 1 + !BETWEEN(t2, t5, t4);
                            t6 = Case6 == 1 ? SUC(t5) : PRED(t5);
                        } else {
                            t6 = t6 == t5->Pred ? t5->Suc : t5->Pred;
                            if ((t5 == t1 && t6 == t2) ||
                                (t5 == t2 && t6 == t1))
                                continue;
                            Case6 += 2;
                        }
                    } else if (BETWEEN(t2, t5, t3)) {
                        Case6 = 4 + X6;
                        t6 = X6 == 1 ? SUC(t5) : PRED(t5);
                        if (t6 == t1)
                            continue;
                    } else {
                        if (X6 == 2)
                            break;
                        Case6 = 7;
                        t6 = PRED(t5);
                        if (t6 == t2)
                            continue;
                    }
                    if (FixedOrCommon(t5, t6))
                        continue;
                    G4 = G3 + C(t5, t6);
                    if ((Case6 <= 2 || Case6 == 5 || Case6 == 6) &&
                        !Forbidden(t6, t1) &&
                        (!c || G4 - c(t6, t1) > 0) &&
                        (*Gain = G4 - C(t6, t1)) > 0) {
                        Make3OptMove(t1, t2, t3, t4, t5, t6, Case6);
                        return 0;
                    }
                    if (Backtracking && !Excludable(t5, t6))
                        continue;
                    Breadth6 = 0;
                    /* Choose (t6,t7) as a candidate edge emanating from t6 */
                    for (Nt6 = t6->CandidateSet; (t7 = Nt6->To); Nt6++) {
                        if (t7 == t6->Pred || t7 == t6->Suc ||
                            (t6 == t2 && t7 == t3) ||
                            (t6 == t3 && t7 == t2) ||
                            ((G5 = G4 - Nt6->Cost) <= 0 &&
                             GainCriterionUsed &&
                             ProblemType != HCP && ProblemType != HPP))
                            continue;
                        if (++Breadth6 > MaxBreadth)
                            break;
                        /* Choose t8 as one of t7's two neighbors on the tour */
                        for (X8 = 1; X8 <= 2; X8++) {
                            if (X8 == 1) {
                                Case8 = Case6;
                                t8 = 0;
                                switch (Case6) {
                                case 1:
                                    t8 = BETWEEN(t2, t7,
                                                 t5) ? SUC(t7) : PRED(t7);
                                    break;
                                case 2:
                                    t8 = BETWEEN(t3, t7,
                                                 t6) ? SUC(t7) : PRED(t7);
                                    break;
                                case 3:
                                    if (BETWEEN(t5, t7, t4))
                                        t8 = SUC(t7);
                                    break;
                                case 4:
                                    if (BETWEEN(t2, t7, t5))
                                        t8 = BETWEEN(t2, t7,
                                                     t4) ? SUC(t7) :
                                            PRED(t7);
                                    break;
                                case 5:
                                    t8 = PRED(t7);
                                    break;
                                case 6:
                                    t8 = BETWEEN(t2, t7,
                                                 t3) ? SUC(t7) : PRED(t7);
                                    break;
                                case 7:
                                    if (BETWEEN(t2, t7, t3))
                                        t8 = SUC(t7);
                                    break;
                                }
                                if (t8 == 0)
                                    break;
                            } else {
                                if (Case6 != 3 && Case6 != 4 && Case6 != 7)
                                    break;
                                t8 = t8 == t7->Pred ? t7->Suc : t7->Pred;
                                Case8 += 8;
                            }
                            if (t8 == t1 ||
                                (t7 == t3 && t8 == t4) ||
                                (t7 == t4 && t8 == t3))
                                continue;
                            if (FixedOrCommon(t7, t8))
                                continue;
                            G6 = G5 + C(t7, t8);
                            if (t8 != t1 &&
                                !Forbidden(t8, t1) &&
                                (!c || G6 - c(t8, t1) > 0) &&
                                (*Gain = G6 - C(t8, t1)) > 0) {
                                Make4OptMove(t1, t2, t3, t4, t5, t6, t7,
                                             t8, Case8);
                                return 0;
                            }
                            if (GainCriterionUsed &&
                                G6 - Precision < t8->Cost)
                                continue;
                            if (!Backtracking || Swaps > 0) {
                                if ((G6 > BestG6 ||
                                     (G6 == BestG6 && !Near(t7, t8) &&
                                      Near(T7, T8))) &&
                                    Swaps < MaxSwaps &&
                                    Excludable(t7, t8) &&
                                    !InInputTour(t7, t8)) {
                                    /* Ignore the move if the gain does 
                                       not vary */
                                    if (RestrictedSearch &&
                                        ProblemType != HCP &&
                                        ProblemType != HPP &&
                                        G2 - t4->Pi == G4 - t6->Pi &&
                                        G4 - t6->Pi == G6 - t8->Pi &&
                                        G3 + t5->Pi == G1 + t3->Pi &&
                                        G5 + t7->Pi == G3 + t5->Pi)
                                        continue;
                                    T3 = t3;
                                    T4 = t4;
                                    T5 = t5;
                                    T6 = t6;
                                    T7 = t7;
                                    T8 = t8;
                                    BestCase8 = Case8;
                                    BestG6 = G6;
                                }
                            } else if (MaxSwaps > 0) {
                                GainType G = G6;
                                Node *t = t8;
                                Make4OptMove(t1, t2, t3, t4, t5, t6, t7,
                                             t8, Case8);
                                Exclude(t1, t2);
                                Exclude(t3, t4);
                                Exclude(t5, t6);
                                Exclude(t7, t8);
                                while ((t =
                                        BestSubsequentMove(t1, t, &G,
                                                           Gain)));
                                if (*Gain > 0)
                                    return 0;
                                RestoreTour();
                                if (t2 != SUC(t1))
                                    Reversed ^= 1;
                            }
                        }
                    }
                }
            }
        }
    }
    *Gain = 0;
    if (T8) {
        /* Make the best 4-opt move */
        Make4OptMove(t1, t2, T3, T4, T5, T6, T7, T8, BestCase8);
        Exclude(t1, t2), Exclude(T3, T4);
        Exclude(T5, T6);
        Exclude(T7, T8);
        *G0 = BestG6;
    }
    return T8;
}

/*
   Below is shown the use of the variables X4, Case6, Case8 and Case10 to 
   discriminate between the 20 cases considered by the algorithm. 

   The notation

        ab-

   is used for a subtour that starts with the edge (ta,tb). For example 
   the tour 

        12-43-

   contains the edges (t1,t2) and (t4,t3), in that order. 

   X4 = 1:
       12-43-
       Case6 = 1: 
           12-56-43-
           Case8 = 1: 
               12-78-56-43-, 12-56-87-43-, 12-56-43-87-
       Case6 = 2:   
           12-43-65-
           Case8 = 2: 
               12-87-43-65-, 12-43-78-65-, 12-43-65-87-
       Case6 = 3: 
           12-65-43-
           Case8 = 3:
               12-65-78-43-
           Case8 = 11:
               12-65-87-43-
       Case6 = 4: 
           12-43-56-
           Case8 = 4: 
               12-78-43-56, 12-43-87-56
           Case8 = 12:
               12-87-43-56-
   X4 = 2:
       12-34-
       Case6 = 5: 
           12-56-34-
           Case8 = 5: 
               12-87-56-34-, 12-56-87-34-, 12-56-34-87-
           Case8 = 13:
               12-56-87-34-
        Case6 = 6: 
           12-65-34-
           Case8 = 6: 
               12-78-65-34-, 12-65-34-87-
           Case8 = 14:
               12-65-87-34-        
       Case6 = 7: 
           12-34-65-
           Case8 = 7: 
               12-78-34-65-
           Case8 = 15:
               12-87-34-65-
*/
//################################### Best4OptMove.c end ###################################

//################################### Best5OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Best5OptMove function makes sequential edge exchanges. If possible, it
 * makes an r-opt move (r <= 5) that improves the tour. Otherwise, it makes
 * the most promising 5-opt move that fulfils the positive gain criterion.
 * To prevent an infinity chain of moves the last edge in a 5-opt move must
 * not previously have been included in the chain.
 *
 * The edge (t1,t2) is the first edge to be exchanged. G0 is a pointer to the
 * accumulated gain.
 *
 * In case a r-opt move is found that improves the tour, the improvement of
 * the cost is made available to the caller through the parameter Gain.
 * If *Gain > 0, an improvement of the current tour has been found. In this
 * case the function returns 0.
 *
 * Otherwise, the best 5-opt move is made, and a pointer to the node that was
 * connected to t1 (in order to close the tour) is returned. The new
 * accumulated gain is made available to the caller through the parameter G0.
 *
 * The function is called from the LinKernighan function.
 */

/* 
   The algorithm splits the set of possible moves up into a number disjoint
   subsets (called "cases"). When t1, t2, ..., t6 has been chosen, Case6 is
   used to discriminate between 8 cases. When t1, t2, ..., t8 has been chosen,
   Case8 is used to discriminate between 16 cases. When t1, t2, ..., t10 has
   been chosen, Case10 is used to discriminate between 52 cases.

   A description of the cases is given after the code.
*/

Node *Best5OptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
{
    Node *t3, *t4, *t5, *t6 = 0, *t7, *t8 = 0, *t9 = 0, *t10 = 0;
    Node *T3 = 0, *T4 = 0, *T5 = 0, *T6 = 0, *T7 = 0, *T8 = 0, *T9 = 0,
        *T10 = 0;
    Candidate *Nt2, *Nt4, *Nt6, *Nt8;
    GainType G1, G2, G3, G4, G5, G6, G7, G8, BestG8 = MINUS_INFINITY;
    int Case6 = 0, Case8 = 0, Case10 = 0, BestCase10 = 0, X4, X6, X8, X10,
        BTW275 = 0, BTW674 = 0, BTW571 = 0, BTW376 = 0, BTW574 = 0,
        BTW671 = 0, BTW471 = 0, BTW673 = 0, BTW573 = 0, BTW273 = 0;
    int Breadth2 = 0, Breadth4, Breadth6, Breadth8;

    if (t2 != SUC(t1))
        Reversed ^= 1;

    /* Determine (T3,T4,T5,T6,T7,T8,T9,T10) = (t3,t4,t5,t6,t7,t8,t9,t10)
       such that

       G8 = *G0 - C(t2,T3) + C(T3,T4) - C(T4,T5) + C(T5,T6) - C(T6,T7) +
       C(T7,T8) - C(T8,T9) + C(T9,T10)

       is maximum (= BestG8), and (T9,T10) has not previously been included.
       If during this process a legal move with *Gain > 0 is found, then
       make the move and exit Best5OptMove immediately. */

    /* Choose (t2,t3) as a candidate edge emanating from t2 */
    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        if (t3 == t2->Pred || t3 == t2->Suc ||
            ((G1 = *G0 - Nt2->Cost) <= 0 && GainCriterionUsed &&
             ProblemType != HCP && ProblemType != HPP))
            continue;
        if (++Breadth2 > MaxBreadth)
            break;
        /* Choose t4 as one of t3's two neighbors on the tour */
        for (X4 = ProblemType == ATSP ? 2 : 1; X4 <= 2; X4++) {
            t4 = X4 == 1 ? PRED(t3) : SUC(t3);
            if (FixedOrCommon(t3, t4))
                continue;
            G2 = G1 + C(t3, t4);
            if (X4 == 1 && !Forbidden(t4, t1) &&
                (!c || G2 - c(t4, t1) > 0) && (*Gain = G2 - C(t4, t1)) > 0)
            {
                Make2OptMove(t1, t2, t3, t4);
                return 0;
            }
            if (Backtracking && !Excludable(t3, t4))
                continue;
            Breadth4 = 0;
            /* Choose (t4,t5) as a candidate edge emanating from t4 */
            for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {
                if (t5 == t4->Pred || t5 == t4->Suc ||
                    ((G3 = G2 - Nt4->Cost) <= 0 && GainCriterionUsed &&
                     ProblemType != HCP && ProblemType != HPP))
                    continue;
                if (++Breadth4 > MaxBreadth)
                    break;
                /* Choose t6 as one of t5's two neighbors on the tour */
                for (X6 = 1; X6 <= 2; X6++) {
                    if (X4 == 1) {
                        if (X6 == 1) {
                            Case6 = 1 + !BETWEEN(t2, t5, t4);
                            t6 = Case6 == 1 ? SUC(t5) : PRED(t5);
                        } else {
                            t6 = t6 == t5->Pred ? t5->Suc : t5->Pred;
                            if ((t5 == t1 && t6 == t2) ||
                                (t5 == t2 && t6 == t1))
                                continue;
                            Case6 += 2;
                        }
                    } else if (BETWEEN(t2, t5, t3)) {
                        Case6 = 4 + X6;
                        t6 = X6 == 1 ? SUC(t5) : PRED(t5);
                        if (t6 == t1)
                            continue;
                    } else {
                        Case6 = 6 + X6;
                        t6 = X6 == 1 ? PRED(t5) : SUC(t5);
                        if (t6 == t2)
                            continue;
                    }
                    if (FixedOrCommon(t5, t6))
                        continue;
                    G4 = G3 + C(t5, t6);
                    if ((Case6 <= 2 || Case6 == 5 || Case6 == 6) &&
                        !Forbidden(t6, t1) &&
                        (!c || G4 - c(t6, t1) > 0) &&
                        (*Gain = G4 - C(t6, t1)) > 0) {
                        Make3OptMove(t1, t2, t3, t4, t5, t6, Case6);
                        return 0;
                    }
                    if (Backtracking && !Excludable(t5, t6))
                        continue;
                    Breadth6 = 0;
                    /* Choose (t6,t7) as a candidate edge emanating from t6 */
                    for (Nt6 = t6->CandidateSet; (t7 = Nt6->To); Nt6++) {
                        if (t7 == t6->Pred || t7 == t6->Suc ||
                            (t6 == t2 && t7 == t3) ||
                            (t6 == t3 && t7 == t2) ||
                            ((G5 = G4 - Nt6->Cost) <= 0 &&
                             GainCriterionUsed &&
                             ProblemType != HCP && ProblemType != HPP))
                            continue;
                        if (++Breadth6 > MaxBreadth)
                            break;
                        /* Choose t8 as one of t7's two neighbors on the tour */
                        for (X8 = 1; X8 <= 2; X8++) {
                            if (X8 == 1) {
                                Case8 = Case6;
                                switch (Case6) {
                                case 1:
                                    if ((BTW275 = BETWEEN(t2, t7, t5)))
                                        t8 = SUC(t7);
                                    else {
                                        t8 = PRED(t7);
                                        BTW674 = BETWEEN(t6, t7, t4);
                                    }
                                    break;
                                case 2:
                                    if ((BTW376 = BETWEEN(t3, t7, t6)))
                                        t8 = SUC(t7);
                                    else {
                                        t8 = PRED(t7);
                                        BTW571 = BETWEEN(t5, t7, t1);
                                    }
                                    break;
                                case 3:
                                    t8 = SUC(t7);
                                    BTW574 = BETWEEN(t5, t7, t4);
                                    break;
                                case 4:
                                    if ((BTW671 = BETWEEN(t6, t7, t1)))
                                        t8 = PRED(t7);
                                    else
                                        t8 = BETWEEN(t2, t7,
                                                     t4) ? SUC(t7) :
                                            PRED(t7);
                                    break;
                                case 5:
                                    t8 = PRED(t7);
                                    BTW471 = BETWEEN(t4, t7, t1);
                                    if (!BTW471)
                                        BTW673 = BETWEEN(t6, t7, t3);
                                    break;
                                case 6:
                                    if ((BTW471 = BETWEEN(t4, t7, t1)))
                                        t8 = PRED(t7);
                                    else {
                                        t8 = SUC(t7);
                                        BTW573 = BETWEEN(t5, t7, t3);
                                    }
                                    break;
                                case 7:
                                case 8:
                                    t8 = SUC(t7);
                                    BTW273 = BETWEEN(t2, t7, t3);
                                    break;
                                }
                            } else {
                                t8 = t8 == t7->Pred ? t7->Suc : t7->Pred;
                                Case8 += 8;
                            }
                            if ((t7 == t1 && t8 == t2) ||
                                (t7 == t2 && t8 == t1) ||
                                (t7 == t3 && t8 == t4) ||
                                (t7 == t4 && t8 == t3))
                                continue;
                            if (FixedOrCommon(t7, t8))
                                continue;
                            if (Case6 == 3 && !BTW574 &&
                                (X8 == 1) == BETWEEN(t3, t7, t1))
                                continue;
                            if (Case6 == 4 && BTW671 && X8 == 2)
                                break;
                            if (Case6 == 7 && !BTW273 &&
                                (X8 == 1) == BETWEEN(t5, t7, t1))
                                continue;
                            if (Case6 == 8 && !BTW273 &&
                                !BETWEEN(t4, t7, t5))
                                break;
                            G6 = G5 + C(t7, t8);
                            if (t8 != t1 &&
                                (Case6 == 3 ? BTW574 :
                                 Case6 == 4 ? !BTW671 :
                                 Case6 == 7 ? BTW273 :
                                 Case6 != 8 && X8 == 1) &&
                                !Forbidden(t8, t1) &&
                                (!c || G6 - c(t8, t1) > 0) &&
                                (*Gain = G6 - C(t8, t1)) > 0) {
                                Make4OptMove(t1, t2, t3, t4, t5, t6, t7,
                                             t8, Case8);
                                return 0;
                            }
                            if (Backtracking && !Excludable(t7, t8))
                                continue;
                            Breadth8 = 0;
                            /* Choose (t8,t9) as a candidate edge emanating 
                               from t8 */
                            for (Nt8 = t8->CandidateSet; (t9 = Nt8->To);
                                 Nt8++) {
                                if (t9 == t8->Pred || t9 == t8->Suc ||
                                    t9 == t1 ||
                                    (t8 == t2 && t9 == t3) ||
                                    (t8 == t3 && t9 == t2) ||
                                    (t8 == t4 && t9 == t5) ||
                                    (t8 == t5 && t9 == t4) ||
                                    ((G7 = G6 - Nt8->Cost) <= 0 &&
                                     GainCriterionUsed &&
                                     ProblemType != HCP &&
                                     ProblemType != HPP))
                                    continue;
                                if (++Breadth8 > MaxBreadth)
                                    break;
                                /* Choose t10 as one of t9's two neighbors
                                   on the tour */
                                for (X10 = 1; X10 <= 2; X10++) {
                                    if (X10 == 1) {
                                        t10 = 0;
                                        switch (Case8) {
                                        case 1:
                                            t10 = (BTW275 ?
                                                   BETWEEN(t8, t9, t5)
                                                   || BETWEEN(t3, t9,
                                                              t1) : BTW674
                                                   ? BETWEEN(t7, t9,
                                                             t1) :
                                                   BETWEEN(t7, t9,
                                                           t5)) ? PRED(t9)
                                                : SUC(t9);
                                            Case10 = 22;
                                            break;
                                        case 2:
                                            t10 = (BTW376 ?
                                                   BETWEEN(t8, t9, t4) :
                                                   BTW571 ?
                                                   BETWEEN(t7, t9, t1)
                                                   || BETWEEN(t3, t9,
                                                              t6) :
                                                   BETWEEN(t7, t9,
                                                           t1)) ? PRED(t9)
                                                : SUC(t9);
                                            Case10 = 23;
                                            break;
                                        case 3:
                                            if (BTW574) {
                                                t10 = BETWEEN(t5, t9, t1) ?
                                                    PRED(t9) : SUC(t9);
                                                Case10 = 24;
                                                break;
                                            }
                                            if (!BETWEEN(t5, t9, t4))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 1;
                                            break;
                                        case 4:
                                            if (BTW671) {
                                                if (!BETWEEN(t2, t9, t5))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 2;
                                                break;
                                            }
                                            t10 = BETWEEN(t6, t9, t4) ?
                                                PRED(t9) : SUC(t9);
                                            Case10 = 25;
                                            break;
                                        case 5:
                                            t10 = (BTW471 ?
                                                   BETWEEN(t7, t9, t1) :
                                                   BTW673 ?
                                                   BETWEEN(t7, t9, t5) :
                                                   BETWEEN(t4, t9, t1)
                                                   || BETWEEN(t7, t9,
                                                              t5)) ?
                                                PRED(t9) : SUC(t9);
                                            Case10 = 26;
                                            break;
                                        case 6:
                                            t10 = (BTW471 ?
                                                   BETWEEN(t7, t9, t3) :
                                                   BTW573 ?
                                                   BETWEEN(t8, t9, t6) :
                                                   BETWEEN(t4, t9, t1)
                                                   || BETWEEN(t8, t9,
                                                              t6)) ?
                                                PRED(t9) : SUC(t9);
                                            Case10 = 27;
                                            break;
                                        case 7:
                                            if (BTW273) {
                                                t10 = BETWEEN(t5, t9, t3) ?
                                                    PRED(t9) : SUC(t9);
                                                Case10 = 28;
                                                break;
                                            }
                                            if (!BETWEEN(t2, t9, t3))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 3;
                                            break;
                                        case 8:
                                            if (BTW273) {
                                                if (!BETWEEN(t4, t9, t5))
                                                    break;
                                                Case10 = 4;
                                            } else {
                                                if (!BETWEEN(t2, t9, t3))
                                                    break;
                                                Case10 = 5;
                                            }
                                            t10 = SUC(t9);
                                            break;
                                        case 9:
                                            if (BTW275) {
                                                if (!BETWEEN(t7, t9, t4))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 6;
                                                break;
                                            }
                                            if (!BTW674) {
                                                if (!BETWEEN(t2, t9, t7))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 7;
                                                break;
                                            }
                                            if (!BETWEEN(t6, t9, t7))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 8;
                                            break;
                                        case 10:
                                            if (BTW376) {
                                                if (!BETWEEN(t7, t9, t6))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 9;
                                                break;
                                            }
                                            if (BTW571) {
                                                if (!BETWEEN(t2, t9, t7))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 10;
                                                break;
                                            }
                                            if (!BETWEEN(t3, t9, t6) &&
                                                !BETWEEN(t2, t9, t7))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 11;
                                            break;
                                        case 11:
                                            if (BTW574) {
                                                t10 = BETWEEN(t3, t9, t1) ?
                                                    PRED(t9) : SUC(t9);
                                                Case10 = 29;
                                                break;
                                            }
                                            if (!BETWEEN(t5, t9, t4))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 12;
                                            break;
                                        case 12:
                                            t10 = BETWEEN(t3, t9, t1) ?
                                                PRED(t9) : SUC(t9);
                                            Case10 = 30;
                                            break;
                                        case 13:
                                            if (BTW471) {
                                                if (!BETWEEN(t2, t9, t7))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 13;
                                                break;
                                            }
                                            if (BTW673) {
                                                if (!BETWEEN(t6, t9, t7))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 14;
                                                break;
                                            }
                                            if (!BETWEEN(t6, t9, t3) &&
                                                !BETWEEN(t2, t9, t7))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 15;
                                            break;
                                        case 14:
                                            if (BTW471) {
                                                if (!BETWEEN(t2, t9, t7))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 16;
                                                break;
                                            }
                                            if (BTW573) {
                                                if (!BETWEEN(t7, t9, t3) &&
                                                    !BETWEEN(t2, t9, t6))
                                                    break;
                                                t10 = SUC(t9);
                                                Case10 = 17;
                                                break;
                                            }
                                            if (!BETWEEN(t7, t9, t6))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 18;
                                            break;
                                        case 15:
                                            if (BTW273) {
                                                t10 = BETWEEN(t5, t9, t1) ?
                                                    PRED(t9) : SUC(t9);
                                                Case10 = 31;
                                                break;
                                            }
                                            if (!BETWEEN(t2, t9, t3))
                                                break;
                                            t10 = SUC(t9);
                                            Case10 = 19;
                                            break;
                                        case 16:
                                            if (BTW273) {
                                                if (!BETWEEN(t4, t9, t5))
                                                    break;
                                                Case10 = 20;
                                            } else {
                                                if (!BETWEEN(t2, t9, t3))
                                                    break;
                                                Case10 = 21;
                                            }
                                            t10 = SUC(t9);
                                            break;
                                        }
                                        if (!t10)
                                            break;
                                    } else {
                                        if (Case10 >= 22)
                                            continue;
                                        Case10 += 31;
                                        t10 =
                                            t10 ==
                                            t9->Pred ? t9->Suc : t9->Pred;
                                    }
                                    if (t10 == t1 ||
                                        (t9 == t3 && t10 == t4) ||
                                        (t9 == t4 && t10 == t3) ||
                                        (t9 == t5 && t10 == t6) ||
                                        (t9 == t6 && t10 == t5))
                                        continue;
                                    if (FixedOrCommon(t9, t10))
                                        continue;
                                    G8 = G7 + C(t9, t10);
                                    if (!Forbidden(t10, t1) &&
                                        (!c || G8 - c(t10, t1) > 0) &&
                                        (*Gain = G8 - C(t10, t1)) > 0) {
                                        Make5OptMove(t1, t2, t3, t4, t5,
                                                     t6, t7, t8, t9, t10,
                                                     Case10);
                                        return 0;
                                    }
                                    if (GainCriterionUsed &&
                                        G8 - Precision < t10->Cost)
                                        continue;
                                    if (!Backtracking || Swaps > 0) {
                                        if ((G8 > BestG8 ||
                                             (G8 == BestG8
                                              && !Near(t9, t10)
                                              && Near(T9, T10)))
                                            && Swaps < MaxSwaps
                                            && Excludable(t9, t10)
                                            && !InInputTour(t9, t10)) {
                                            /* Ignore the move if the gain
                                               does not vary */
                                            if (RestrictedSearch &&
                                                ProblemType != HCP &&
                                                ProblemType != HPP &&
                                                G2 - t4->Pi == G4 - t6->Pi
                                                && G4 - t6->Pi ==
                                                G6 - t8->Pi
                                                && G6 - t8->Pi ==
                                                G8 - t10->Pi
                                                && G3 + t5->Pi ==
                                                G1 + t3->Pi
                                                && G5 + t7->Pi ==
                                                G3 + t5->Pi
                                                && G7 + t9->Pi ==
                                                G5 + t7->Pi)
                                                continue;
                                            T3 = t3;
                                            T4 = t4;
                                            T5 = t5;
                                            T6 = t6;
                                            T7 = t7;
                                            T8 = t8;
                                            T9 = t9;
                                            T10 = t10;
                                            BestCase10 = Case10;
                                            BestG8 = G8;
                                        }
                                    } else if (MaxSwaps > 0) {
                                        GainType G = G8;
                                        Node *t = t10;
                                        Make5OptMove(t1, t2, t3, t4, t5,
                                                     t6, t7, t8, t9, t10,
                                                     Case10);
                                        Exclude(t1, t2);
                                        Exclude(t3, t4);
                                        Exclude(t5, t6);
                                        Exclude(t7, t8);
                                        Exclude(t9, t10);
                                        while ((t =
                                                BestSubsequentMove(t1, t,
                                                                   &G,
                                                                   Gain)));
                                        if (*Gain > 0)
                                            return 0;
                                        RestoreTour();
                                        if (t2 != SUC(t1))
                                            Reversed ^= 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    *Gain = 0;
    if (T10) {
        /* Make the best 5-opt move */
        Make5OptMove(t1, t2, T3, T4, T5, T6, T7, T8, T9, T10, BestCase10);
        Exclude(t1, t2);
        Exclude(T3, T4);
        Exclude(T5, T6);
        Exclude(T7, T8);
        Exclude(T9, T10);
        *G0 = BestG8;
    }
    return T10;
}

/*
   Below is shown the use of the variables X4, Case6, Case8 and Case10 to
   discriminate between the 148 cases considered by the algorithm.

   The notation

        ab-

   is used for a subtour that starts with the edge (ta,tb). For example
   the tour

        12-43-

   contains the edges (t1,t2) and (t4,t3), in that order.

   X4 = 1:
       12-43-
       Case6 = 1:
           12-56-43-
           Case8 = 1:
               12-78-56-43-, 12-56-87-43-, 12-56-43-87-
               Case10 = 22:
                   12-910-78-56-43-, 12-78-109-56-43-, 12-78-56-910-43-, 12-78-56-43-109-
                   12-910-56-87-43-, 12-56-910-87-43-, 12-56-87-109-43-, 12-56-87-43-109-
                   12-109-56-43-87-, 12-56-910-43-87-, 12-56-43-910-87-, 12-56-43-87-109-
           Case8 = 9:
               12-87-56-43-, 12-56-78-43-, 12-56-43-78-
               Case10 = 6:
                   12-87-910-56-43-, 12-87-56-910-43-
               Case10 = 7:
                   12-910-56-43-78-, 12-56-910-43-78-, 12-56-43-910-78-
               Case10 = 8:
                   12-56-910-78-43-
               Case10 = 37:
                   12-87-109-56-43-, 12-87-56-109-43-
               Case10 = 38:
                   12-109-56-43-78-, 12-56-109-43-78-, 12-56-43-109-78-
               Case10 = 39:
                   12-56-109-78-43-
       Case6 = 2:
           12-43-65-
           Case8 = 2:
               12-87-43-65-, 12-43-78-65-, 12-43-65-87-
               Case10 = 23:
                   12-910-87-43-65-, 12-87-910-43-65-, 12-87-43-910-65-, 12-87-43-65-910-
                   12-109-43-78-65-, 12-43-910-78-65-, 12-43-78-109-65-, 12-43-78-65-109-
                   12-910-43-65-87-, 12-43-109-65-87-, 12-43-65-910-87-, 12-43-65-87-109-
           Case8 = 10:
               12-78-43-65-, 12-43-87-65-, 12-43-65-78-
               Case10 = 9:
                   12-43-87-910-65-
               Case10 = 10:
                   12-910-43-65-78-, 12-43-910-65-78-, 12-43-65-910-78-
               Case10 = 11:
                   12-910-78-43-65-, 12-78-43-910-43-
               Case10 = 40:
                   12-43-87-109-65-
               Case10 = 41:
                   12-109-43-65-78-, 12-43-109-65-78-, 12-43-65-109-78-
               Case10 = 42:
                   12-109-78-43-65-, 12-78-43-109-43-
       Case6 = 3:
           12-65-43-
           Case8 = 3:
               12-78-65-43-, 12-65-78-43-
               Case10 = 24:
                   12-910-65-78-43-, 12-65-109-78-43-, 12-65-78-109-43-, 12-65-78-43-109-
               Case10 = 1:
                   12-78-65-910-43-
               Case10 = 32:
                   12-78-65-109-43
           Case8 = 11:
               12-65-87-43-, 12-65-43-87-
               Case10 = 29:
                   12-910-65-87-43-, 12-65-910-87-43-, 12-65-87-910-43-, 12-65-87-43-109-
               Case10 = 12:
                   12-65-910-43-87-
               Case10 = 43:
                   12-65-109-43-87-
       Case6 = 4:
           12-43-56-
           Case8 = 4:
               12-78-43-56, 12-43-87-56, 12-43-56-87-
               Case10 = 2:
                   12-910-43-56-87-, 12-43-910-56-87-
               Case10 = 25:
                   12-109-78-43-56-, 12-78-109-43-56-, 12-67-43-910-56-, 12-67-43-56-109-
                   12-109-43-87-56-, 12-43-910-87-56-, 12-43-87-910-56-, 12-43-87-56-109-
               Case10 = 33:
                   12-109-43-56-87-, 12-43-109-56-87-
           Case8 = 12:
               12-87-43-56-, 12-43-78-56-
           Case10 = 30:
               12-910-87-43-56-, 12-87-910-43-56-, 12-87-43-109-56-, 12-87-43-56-109-
               12-910-43-78-56-, 12-43-109-78-56-, 12-43-78-109-56-, 12-43-78-56-109-
   X4 = 2:
       12-34-
       Case6 = 5:
           12-56-34-
           Case8 = 5: 
               12-87-56-34-, 12-56-87-34-, 12-56-34-87-
               Case10 = 26:
                   12-910-87-56-34-, 12-87-109-56-34-, 12-87-56-910-34-, 12-87-56-34-109-
                   12-109-56-87-34-, 12-56-910-87-34-, 12-56-87-109-34-, 12-56-87-34-109-
                   12-910-56-34-87-, 12-56-910-34-87-, 12-56-34-910-87-, 12-56-34-87-109-
           Case8 = 13:
               12-78-56-34-, 12-56-78-34-, 12-56-34-78-
               Case10 = 13:
                   12-910-56-34-78-, 12-56-910-34-78-, 12-56-34-910-78-
               Case10 = 14:
                   12-56-910-78-34-
               Case10 = 15:
                   12-910-78-56-34-, 12-78-56-910-34-
               Case10 = 44:
                   12-109-56-34-78-, 12-56-109-34-78-, 12-56-34-109-78-
               Case10 = 45:
                   12-56-109-78-34-
               Case10 = 46:
                   12-109-78-56-34-, 12-78-56-109-34-
       Case6 = 6:
           12-65-34-
           Case8 = 6: 
               12-87-65-34-, 12-65-78-34-, 12-65-34-87-
               Case10 = 27:
                   12-910-87-65-34-, 12-87-109-65-34-, 12-87-65-910-34-, 12-87-65-34-109-
                   12-109-65-78-34-, 12-65-910-78-34-, 12-65-78-109-34-, 12-65-78-34-109-
                   12-910-65-34-87-, 12-65-910-34-87-, 12-65-34-910-87-, 12-65-34-87-109-
           Case8 = 14:
               12-87-65-34-, 12-65-87-34-, 12-65-34-78-
               Case10 = 16:
                   12-910-65-34-78-, 12-65-910-34-78-, 12-65-34-910-78-
               Case10 = 17:
                   12-910-65-87-34-, 12-65-87-910-34-
               Case10 = 18:
                   12-87-910-65-34-
               Case10 = 47:
                  12-109-65-34-78-, 12-65-109-34-78-, 12-65-34-109-78-
               Case10 = 48:
                  12-109-65-87-34-, 12-65-87-109-34-
               Case10 = 49:
                  12-87-109-65-34-
       Case6 = 7:
           12-34-65-
            Case8 = 7:
                12-78-34-65-, 12-34-78-65-
                Case10 = 28:
                    12-109-78-34-65-, 12-78-109-34-65-, 12-78-34-910-34-, 12-78-34-65-109-
                Case10 = 3:
                    12-910-34-78-65-
                Case10 = 34:
                    12-109-34-78-65-
            Case8 = 15:
                12-87-34-65-, 12-34-87-65-
                Case10 = 31:
                    12-910-87-34-65-, 12-87-910-34-65-, 12-87-34-910-65-, 12-87-34-65-109-
                Case10 = 19:
                    12-910-34-87-65-
                Case10 = 50:
                    12-109-34-87-65-
       Case6 = 8:
           12-34-56-
           Case8 = 8: 
               12-78-34-56-, 12-34-78-56-
               Case10 = 4:
                   12-78-34-910-56-
               Case10 = 5:
                   12-910-34-78-56-
               Case10 = 35:
                   12-78-34-109-56-
               Case10 = 36:
                   12-109-34-78-56-
           Case8 = 16:
               12-87-34-56-, 12-34-87-56-
               Case10 = 20:
                   12-87-34-910-56-
               Case10 = 21:
                   12-910-34-87-56-
               Case10 = 51:
                   12-87-34-109-56-
               Case10 = 52:
                   12-109-34-87-56-
*/
//################################### Best5OptMove.c end ###################################

//################################### BestKOptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"
#include "Sequence.h"

/*
 * The BestKOptMove function makes edge exchanges. If possible, it makes a 
 * r-opt move (r >= 2) that improves the tour. Otherwise, it makes the most 
 * promising sequential K-opt move that fulfils the positive gain criterion. 
 * To prevent an infinity chain of moves the last edge in a K-opt move must
 * not previously have been included in the chain. 
 *
 * The edge (t[1],t[2]) is the first edge to be exchanged. G0 is a pointer to 
 * the accumulated gain.
 *
 * In case a K-opt move is found that improves the tour, the improvement of 
 * the cost is made available to the caller through the parameter Gain. 
 * If *Gain > 0, an improvement of the current tour has been found. In this
 * case the function returns 0.
 *
 * Otherwise, the best K-opt move is made, and a pointer to the node that was 
 * connected to t[1] (in order to close the tour) is returned. The new 
 * accumulated gain is made available to the caller through the parameter G0. 
 *
 * The function is called from the LinKernighan function. 
 */

static GainType BestG2;

static GainType BestKOptMoveRec(int k, GainType G0);

Node *BestKOptMove(Node * t1, Node * t2, GainType * G0, GainType * Gain)
{
    K = Swaps == 0 ? MoveType : SubsequentMoveType;
    *Gain = 0;
    t[1] = t1;
    t[2] = t2;
    T[2 * K] = 0;
    BestG2 = MINUS_INFINITY;

    /* 
     * Determine (T[3],T[4], ..., T[2K]) = (t[3],t[4], ..., t[2K])
     * such that
     *
     *     G[2 * K] = *G0 - C(t[2],T[3]) + C(T[3],T[4])
     *                    - C(T[4],T[5]) + C(T[5],T[6])
     *                      ...
     *                    - C(T[2K-3],T[2K-2]) + C(T[2K-1],T[2K])
     *
     * is maximum, and (T[2K-1],T[2K]) has not previously been included.
     * If during this process a legal move with *Gain > 0 is found, then 
     * make the move and exit BestKOptMove immediately.
     */

    MarkDeleted(t1, t2);
    *Gain = BestKOptMoveRec(2, *G0);
    UnmarkDeleted(t1, t2);

    if (*Gain <= 0 && T[2 * K]) {
        int i;
        memcpy(t + 1, T + 1, 2 * K * sizeof(Node *));
        for (i = 2; i < 2 * K; i += 2)
            incl[incl[i] = i + 1] = i;
        incl[incl[1] = 2 * K] = 1;
        MakeKOptMove(K);
        for (i = 1; i < 2 * K; i += 2)
            Exclude(T[i], T[i + 1]);
        *G0 = BestG2;
        return T[2 * K];
    }
    return 0;
}

static GainType BestKOptMoveRec(int k, GainType G0)
{
    Candidate *Nt2;
    Node *t1, *t2, *t3, *t4;
    GainType G1, G2, G3, Gain;
    int X4, i;
    int Breadth2 = 0;

    t1 = t[1];
    t2 = t[i = 2 * k - 2];
    incl[incl[i] = i + 1] = i;
    incl[incl[1] = i + 2] = 1;
    /* Choose (t2,t3) as a candidate edge emanating from t2 */
    for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
        if (t3 == t2->Pred || t3 == t2->Suc ||
            ((G1 = G0 - Nt2->Cost) <= 0 && GainCriterionUsed &&
             ProblemType != HCP && ProblemType != HPP) ||
            Added(t2, t3))
            continue;
        if (++Breadth2 > MaxBreadth)
            break;
        MarkAdded(t2, t3);
        t[2 * k - 1] = t3;
        G[2 * k - 2] = G1 + t3->Pi;
        /* Choose t4 as one of t3's two neighbors on the tour */
        for (X4 = 1; X4 <= 2; X4++) {
            t4 = X4 == 1 ? PRED(t3) : SUC(t3);
            if (FixedOrCommon(t3, t4) || Deleted(t3, t4))
                continue;
            t[2 * k] = t4;
            G2 = G1 + C(t3, t4);
            G3 = MINUS_INFINITY;
            if (t4 != t1 && !Forbidden(t4, t1) && !Added(t4, t1) &&
                (!c || G2 - c(t4, t1) > 0) &&
                (G3 = G2 - C(t4, t1)) > 0 && FeasibleKOptMove(k)) {
                UnmarkAdded(t2, t3);
                MakeKOptMove(k);
                return G3;
            }
            if (Backtracking && !Excludable(t3, t4))
                continue;
            MarkDeleted(t3, t4);
            G[2 * k - 1] = G2 - t4->Pi;
            if (k < K) {
                if ((Gain = BestKOptMoveRec(k + 1, G2)) > 0) {
                    UnmarkAdded(t2, t3);
                    UnmarkDeleted(t3, t4);
                    return Gain;
                }
                incl[incl[1] = 2 * k] = 1;
            }
            if (t4 != t1 && !Forbidden(t4, t1) &&
                k + 1 < NonsequentialMoveType &&
                PatchingC >= 2 && PatchingA >= 1 &&
                (Swaps == 0 || SubsequentPatching)) {
                if (G3 == MINUS_INFINITY)
                    G3 = G2 - C(t4, t1);
                if ((PatchingCRestricted ? G3 > 0 && IsCandidate(t4, t1) :
                     PatchingCExtended ? G3 > 0
                     || IsCandidate(t4, t1) : G3 > 0)
                    && (Gain = PatchCycles(k, G3)) > 0) {
                    UnmarkAdded(t2, t3);
                    UnmarkDeleted(t3, t4);
                    return Gain;
                }
            }
            UnmarkDeleted(t3, t4);
            if (k == K && t4 != t1 && t3 != t1 && G3 <= 0 &&
                !Added(t4, t1) &&
                (!GainCriterionUsed || G2 - Precision >= t4->Cost)) {
                if (!Backtracking || Swaps > 0) {
                    if ((G2 > BestG2 ||
                         (G2 == BestG2 && !Near(t3, t4) &&
                          Near(T[2 * K - 1], T[2 * K]))) &&
                        Swaps < MaxSwaps &&
                        Excludable(t3, t4) && !InInputTour(t3, t4)) {
                        if (RestrictedSearch && K > 2 &&
                            ProblemType != HCP && ProblemType != HPP) {
                            /* Ignore the move if the gain does not vary */
                            G[0] = G[2 * K - 2];
                            G[1] = G[2 * K - 1];
                            for (i = 2 * K - 3; i >= 2; i--)
                                if (G[i] != G[i % 2])
                                    break;
                            if (i < 2)
                                continue;
                        }
                        if (FeasibleKOptMove(K)) {
                            BestG2 = G2;
                            memcpy(T + 1, t + 1, 2 * K * sizeof(Node *));
                        }
                    }
                } else if (MaxSwaps > 0 && FeasibleKOptMove(K)) {
                    Node *SUCt1 = SUC(t1);
                    MakeKOptMove(K);
                    for (i = 1; i < 2 * k; i += 2) {
                        Exclude(t[i], t[i + 1]);
                        UnmarkDeleted(t[i], t[i + 1]);
                    }
                    for (i = 2; i < 2 * k; i += 2)
                        UnmarkAdded(t[i], t[i + 1]);
                    memcpy(tSaved + 1, t + 1, 2 * k * sizeof(Node *));
                    while ((t4 = BestSubsequentMove(t1, t4, &G2, &Gain)));
                    if (Gain > 0) {
                        UnmarkAdded(t2, t3);
                        return Gain;
                    }
                    RestoreTour();
                    K = k;
                    memcpy(t + 1, tSaved + 1, 2 * K * sizeof(Node *));
                    for (i = 1; i < 2 * K - 2; i += 2)
                        MarkDeleted(t[i], t[i + 1]);
                    for (i = 2; i < 2 * K; i += 2)
                        MarkAdded(t[i], t[i + 1]);
                    for (i = 2; i < 2 * K; i += 2)
                        incl[incl[i] = i + 1] = i;
                    incl[incl[1] = 2 * K] = 1;
                    if (SUCt1 != SUC(t1))
                        Reversed ^= 1;
                    T[2 * K] = 0;
                }
            }
        }
        UnmarkAdded(t2, t3);
        if (t3 == t1)
            continue;
        /* Try to delete an added edge, (_,t3) or (t3,_) */
        for (i = 2 * k - 4; i >= 2; i--) {
            if (t3 == t[i]) {
                t4 = t[i ^ 1];
                if (t4 == t1 || Forbidden(t4, t1) || FixedOrCommon(t3, t4)
                    || Added(t4, t1))
                    continue;
                G2 = G1 + C(t3, t4);
                if ((!c || G2 - c(t4, t1) > 0)
                    && (Gain = G2 - C(t4, t1)) > 0) {
                    incl[incl[i ^ 1] = 1] = i ^ 1;
                    incl[incl[i] = 2 * k - 2] = i;
                    if (FeasibleKOptMove(k - 1)) {
                        MakeKOptMove(k - 1);
                        return Gain;
                    }
                    incl[incl[i ^ 1] = i] = i ^ 1;
                }
            }
        }
        incl[1] = 2 * k;
        incl[2 * k - 2] = 2 * k - 1;
    }
    return 0;
}
//################################### BestKOptMove.c end ###################################

//################################### Between.c begin ###################################
#include "LKH.h"

/* 
 * The Between function is used to determine whether a node
 * is between two other nodes with respect to the current
 * orientation. The function is only used if the doubly linked
 * list representation is used for a tour; if the two-la
 * tree representation is used, the function Between_SL is used
 * instead.
 *	
 * Between(ta,tb,tc) returns 1 if node tb is between node ta and tc.
 * Otherwise, 0 is returned.
 *	
 * The function is called from the functions BestMove, Gain23,
 * BridgeGain, Make4OptMove, Make5OptMove and FindPermutation.
 */

int Between(const Node * ta, const Node * tb, const Node * tc)
{
    int a, b = tb->Rank, c;

    if (!Reversed) {
        a = ta->Rank;
        c = tc->Rank;
    } else {
        a = tc->Rank;
        c = ta->Rank;
    }
    return a <= c ? b >= a && b <= c : b >= a || b <= c;
}
//################################### Between.c end ###################################

//################################### Between_SL.c begin ###################################
#include "LKH.h"

/*
 * The Between_SL function is used to determine whether a node is 
 * between two other nodes with respect to the current orientation. 
 * The function is only used if the two-level tree representation 
 * is used for a tour; if the doubly linked list representation is 
 * used, the function Between is used instead.
 * 	
 * Between_SL(a,b,c) returns 1 if node b is between node a and c.
 * Otherwise, 0 is returned.
 * 	
 * The function is called from the functions BestMove, Gain23,
 * BridgeGain, Make4OptMove, Make5OptMove and FindPermutation.
 */

int Between_SL(const Node * ta, const Node * tb, const Node * tc)
{
    const Segment *Pa, *Pb, *Pc;

    if (tb == ta || tb == tc)
        return 1;
    if (ta == tc)
        return 0;
    Pa = ta->Parent;
    Pb = tb->Parent;
    Pc = tc->Parent;
    if (Pa == Pc) {
        if (Pb == Pa)
            return (Reversed == Pa->Reversed) ==
                (ta->Rank < tc->Rank ?
                 tb->Rank > ta->Rank && tb->Rank < tc->Rank :
                 tb->Rank > ta->Rank || tb->Rank < tc->Rank);
        return (Reversed == Pa->Reversed) == (ta->Rank > tc->Rank);
    }
    if (Pb == Pc)
        return (Reversed == Pb->Reversed) == (tb->Rank < tc->Rank);
    if (Pa == Pb)
        return (Reversed == Pa->Reversed) == (ta->Rank < tb->Rank);
    return Reversed !=
        (Pa->Rank < Pc->Rank ?
         Pb->Rank > Pa->Rank && Pb->Rank < Pc->Rank :
         Pb->Rank > Pa->Rank || Pb->Rank < Pc->Rank);
}
//################################### Between_SL.c end ###################################

//################################### Between_SSL.c begin ###################################
#include "LKH.h"

/*
   The Between_SSL function is used to determine whether a segment is 
   between two other segments with respect to the current orientation. 
   The function is only used if the three-level tree representation 
   is used for a tour.
   	
   Between_SSL(a,b,c) returns 1 if segment b is between segment a and c.
   Otherwise, 0 is returned.
   	
   The function is called from the functions BestMove, Gain23,
   BridgeGain, Make4OptMove and Make5OptMove.
*/

int Between_SSL(const Node * ta, const Node * tb, const Node * tc)
{
    const Segment *Pa, *Pb, *Pc;
    const SSegment *PPa, *PPb, *PPc;

    if (tb == ta || tb == tc)
        return 1;
    if (ta == tc)
        return 0;
    Pa = ta->Parent;
    Pb = tb->Parent;
    Pc = tc->Parent;
    PPa = Pa->Parent;
    PPb = Pb->Parent;
    PPc = Pc->Parent;
    if (Pa == Pc) {
        if (Pb == Pa)
            return (Reversed == (Pa->Reversed != PPa->Reversed)) ==
                (ta->Rank < tc->Rank ?
                 tb->Rank > ta->Rank && tb->Rank < tc->Rank :
                 tb->Rank > ta->Rank || tb->Rank < tc->Rank);
        return (Reversed == (Pa->Reversed != PPa->Reversed)) == (ta->Rank >
                                                                 tc->Rank);
    }
    if (Pb == Pc)
        return (Reversed == (Pb->Reversed != PPb->Reversed)) == (tb->Rank <
                                                                 tc->Rank);
    if (Pa == Pb)
        return (Reversed == (Pa->Reversed != PPa->Reversed)) == (ta->Rank <
                                                                 tb->Rank);
    if (PPa == PPc) {
        if (PPb == PPa)
            return (Reversed == PPa->Reversed) ==
                (Pa->Rank < Pc->Rank ?
                 Pb->Rank > Pa->Rank && Pb->Rank < Pc->Rank :
                 Pb->Rank > Pa->Rank || Pb->Rank < Pc->Rank);
        return (Reversed == PPa->Reversed) == (Pa->Rank > Pc->Rank);
    }
    if (PPb == PPc)
        return (Reversed == PPb->Reversed) == (Pb->Rank < Pc->Rank);
    if (PPa == PPb)
        return (Reversed == PPa->Reversed) == (Pa->Rank < Pb->Rank);
    return Reversed !=
        (PPa->Rank < PPc->Rank ?
         PPb->Rank > PPa->Rank && PPb->Rank < PPc->Rank :
         PPb->Rank > PPa->Rank || PPb->Rank < PPc->Rank);
}
//################################### Between_SSL.c end ###################################

//################################### BridgeGain.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The BridgeGain function attempts to improve the tour by making a 
 * non-sequential move. 
 * 
 * The function is called by the Gain23 function. 
 *	
 * For any nonfeasible 2-opt move that would cause the current tour to be 
 * split into two separate tours, BridgGain may be called in order to find 
 * a (nonfeasible) 2- or 3-opt move that reconnects the two separate tours 
 * into a tour which is shorter than the original one. In some cases, the 
 * second move may even be 4-opt.
 * 
 * For any nonfeasible 3-opt move that would cause the current tour to be 
 * split into two separate tours, BridgGain may be called in order to find 
 * a (nonfeasible) 2-opt move that reconnects the two separate tours into 
 * a tour which is shorter than the original one.
 *	
 * The parameters s1, s2, ..., s8 denote the end nodes of edges that are 
 * part of the nonfeasible move suggested by Gain23. The parameter Case6 
 * is used to specify the move type (Case6 = 0 when 2-opt, Case6 = 3, 4 
 * or 7 when 3- or 4-opt). The parameter G  contains the gain achieved by 
 * making the move.
 *	
 * If the composite move results in a shorter tour, then the move is made, 
 * and the function returns the gain achieved.	
 */

GainType
BridgeGain(Node * s1, Node * s2, Node * s3, Node * s4,
           Node * s5, Node * s6, Node * s7, Node * s8, int Case6,
           GainType G)
{
    Node *t1, *t2, *t3, *t4, *t5, *t6, *t7, *t8, *u2 = 0, *u3 = 0;
    Candidate *Nt2, *Nt4, *Nt6;
    GainType G0, G1, G2, G3, G4, G5, G6, Gain;
    int X4;
    int Breadth2, Breadth4, Breadth6;

    /* From the original tour select a segment (u2 --> u3) which contains 
       as few nodes as possible */
    switch (Case6) {
    case 3:
        if (2 * SegmentSize(s5, s4) <= Dimension) {
            u2 = s5;
            u3 = s4;
        } else {
            u2 = s3;
            u3 = s6;
        }
        break;
    case 4:
        if (2 * SegmentSize(s2, s5) <= Dimension) {
            u2 = s2;
            u3 = s5;
        } else {
            u2 = s6;
            u3 = s1;
        }
        break;
    case 0:
    case 7:
        if (2 * SegmentSize(s2, s3) <= Dimension) {
            u2 = s2;
            u3 = s3;
        } else {
            u2 = s4;
            u3 = s1;
        }
    }

    /* Choose t1 between u2 and u3 */
    for (t1 = u2; t1 != u3; t1 = t2) {
        /* Choose t2 as the successor of t1 */
        t2 = SUC(t1);
        if ((t1 == s1 && t2 == s2) ||
            (t1 == s2 && t2 == s1) ||
            (t1 == s3 && t2 == s4) ||
            (t1 == s4 && t2 == s3) ||
            (t1 == s5 && t2 == s6) ||
            (t1 == s6 && t2 == s5) ||
            (t1 == s7 && t2 == s8) ||
            (t1 == s8 && t2 == s7) || FixedOrCommon(t1, t2))
            continue;
        G0 = G + C(t1, t2);
        /* Choose (t2,t3) as a candidate edge emanating from t2. 
           t3 must not be between u2 and u3 */
        Breadth2 = 0;
        for (Nt2 = t2->CandidateSet; (t3 = Nt2->To); Nt2++) {
            if (t3 == t2->Pred || t3 == t2->Suc || BETWEEN(u2, t3, u3))
                continue;
            G1 = G0 - Nt2->Cost;
            if (++Breadth2 > MaxBreadth)
                break;
            /* Choose t4 as one of t3's two neighbors on the tour */
            for (X4 = 1; X4 <= 2; X4++) {
                t4 = X4 == 1 ? SUC(t3) : PRED(t3);
                if (t4 == t2 ||
                    (t3 == s1 && t4 == s2) ||
                    (t3 == s2 && t4 == s1) ||
                    (t3 == s3 && t4 == s4) ||
                    (t3 == s4 && t4 == s3) ||
                    (t3 == s5 && t4 == s6) ||
                    (t3 == s6 && t4 == s5) ||
                    (t3 == s7 && t4 == s8) ||
                    (t3 == s8 && t4 == s7) || FixedOrCommon(t3, t4))
                    continue;
                G2 = G1 + C(t3, t4);
                /* Test if an improvement can be obtained */
                if (!Forbidden(t4, t1) &&
                    (!c || G2 - c(t4, t1) > 0) &&
                    (t4 != s2 || t1 != s3) &&
                    (t4 != s3 || t1 != s2) &&
                    (t4 != s4 || t1 != s5) &&
                    (t4 != s5 || t1 != s4) &&
                    (t4 != s6 || t1 != s7) &&
                    (t4 != s7 || t1 != s6) &&
                    (t4 != s8 || t1 != s1) &&
                    (t4 != s1 || t1 != s8) &&
                    (s8 ||
                     ((t4 != s6 || t1 != s1) &&
                      (t4 != s1 || t1 != s6))) &&
                    (Gain = G2 - C(t4, t1)) > 0) {
                    switch (Case6) {
                    case 0:
                        if (X4 == 1)
                            Swap3(s1, s2, s4, t3, t4, t1, s1, s3, s2);
                        else
                            Swap2(t1, t2, t3, s1, s2, s3);
                        return Gain;
                    case 3:
                        if ((X4 == 1) ==
                            (!BETWEEN(s2, t1, s6) && !BETWEEN(s2, t3, s6)))
                            Swap3(s1, s2, s3, t1, t2, t3, s5, s6, s1);
                        else
                            Swap4(s1, s2, s3, t1, t2, t4, s5, s6, s1, t2,
                                  t4, t1);
                        if (s8)
                            Swap1(s7, s8, s1);
                        return Gain;
                    case 4:
                        if ((X4 == 1) ==
                            (!BETWEEN(s3, t1, s5) && !BETWEEN(s3, t3, s5)))
                            Swap3(s1, s2, s3, t1, t2, t3, s5, s6, s1);
                        else
                            Swap4(s1, s2, s3, t1, t2, t4, s5, s6, s1, t2,
                                  t4, t1);
                        if (s8)
                            Swap1(s7, s8, s1);
                        return Gain;
                    case 7:
                        if ((X4 == 1) ==
                            (!BETWEEN(s4, t1, s6) && !BETWEEN(s4, t3, s6)))
                            Swap3(s5, s6, s1, t1, t2, t3, s3, s4, s5);
                        else
                            Swap4(s5, s6, s1, t1, t2, t4, s3, s4, s5, t2,
                                  t4, t1);
                        if (s8)
                            Swap1(s7, s8, s1);
                        return Gain;
                    }
                }
                /* If BridgeGain has been called with a nonfeasible 2-opt move,
                   then try to find a 3-opt or 4-opt move which, when composed 
                   with the 2-opt move, results in an improvement of the tour */
                if (Case6 != 0)
                    continue;
                Breadth4 = 0;
                /* Choose (t4,t5) as a candidate edge emanating from t4 */
                for (Nt4 = t4->CandidateSet; (t5 = Nt4->To); Nt4++) {
                    if (t5 == t4->Pred || t5 == t4->Suc || t5 == t1
                        || t5 == t2)
                        continue;
                    /* Choose t6 as one of t5's two neighbors on the tour.
                       Only one choice! */
                    t6 = X4 == 1
                        || BETWEEN(u2, t5, u3) ? PRED(t5) : SUC(t5);
                    if ((t5 == s1 && t6 == s2) || (t5 == s2 && t6 == s1)
                        || (t5 == s3 && t6 == s4) || (t5 == s4 && t6 == s3)
                        || FixedOrCommon(t5, t6))
                        continue;
                    G3 = G2 - Nt4->Cost;
                    G4 = G3 + C(t5, t6);
                    if (!Forbidden(t6, t1) &&
                        (!c || G4 - c(t6, t1) > 0) &&
                        (Gain = G4 - C(t6, t1)) > 0) {
                        if (X4 == 1)
                            Swap4(s1, s2, s4, t3, t4, t1, s1, s3, s2, t5,
                                  t6, t1);
                        else
                            Swap3(t1, t2, t3, s1, s2, s3, t5, t6, t1);
                        return Gain;
                    }
                    if (++Breadth4 > MaxBreadth)
                        break;
                    Breadth6 = 0;
                    /* Choose (t7,t8) as a candidate edge emanating from t7.
                       Only one choice! */
                    for (Nt6 = t6->CandidateSet; (t7 = Nt6->To); Nt6++) {
                        if (t7 == t6->Pred || t7 == t6->Suc)
                            continue;
                        /* Choose t8 as one of t7's two neighbors on the tour.
                           Only one choice! */
                        if (X4 == 1)
                            t8 = (BETWEEN(u2, t5, t1) ? BETWEEN(t5, t7, t1)
                                  : BETWEEN(t2, t5, u3) ? BETWEEN(u2, t7,
                                                                  t1)
                                  || BETWEEN(t5, t7, u3) : BETWEEN(SUC(u3),
                                                                   t5,
                                                                   t3) ?
                                  BETWEEN(u2, t7, u3)
                                  || BETWEEN(t5, t7, t3) : !BETWEEN(t4, t7,
                                                                    t6)) ?
                                PRED(t7) : SUC(t7);
                        else
                            t8 = (BETWEEN(u2, t5, t1) ?
                                  !BETWEEN(u2, t7, t6)
                                  && !BETWEEN(t2, t7, u3) : BETWEEN(t2, t5,
                                                                    u3) ?
                                  !BETWEEN(t2, t7, t6) : BETWEEN(SUC(u3),
                                                                 t5,
                                                                 t4) ?
                                  !BETWEEN(SUC(u3), t7, t5)
                                  && !BETWEEN(t3, t7,
                                              PRED(u2)) : !BETWEEN(t3, t7,
                                                                   t5)) ?
                                PRED(t7) : SUC(t7);
                        if (t8 == t1
                            || (t7 == t1 && t8 == t2) || (t7 == t3
                                                          && t8 == t4)
                            || (t7 == t4 && t8 == t3) || (t7 == s1
                                                          && t8 == s2)
                            || (t7 == s2 && t8 == s1) || (t7 == s3
                                                          && t8 == s4)
                            || (t7 == s4 && t8 == s3))
                            continue;
                        if (FixedOrCommon(t7, t8) || Forbidden(t8, t1))
                            continue;
                        G5 = G4 - Nt6->Cost;
                        G6 = G5 + C(t7, t8);
                        /* Test if an improvement can be achieved */
                        if ((!c || G6 - c(t8, t1) > 0) &&
                            (Gain = G6 - C(t8, t1)) > 0) {
                            if (X4 == 1)
                                Swap4(s1, s2, s4, t3, t4, t1, s1, s3, s2,
                                      t5, t6, t1);
                            else
                                Swap3(t1, t2, t3, s1, s2, s3, t5, t6, t1);
                            Swap1(t7, t8, t1);
                            return Gain;
                        }
                        if (++Breadth6 > MaxBreadth)
                            break;
                    }
                }
            }
        }
    }
    /* No improvement has been found */
    return 0;
}
//################################### BridgeGain.c end ###################################

//################################### C.c begin ###################################
#include "LKH.h"

/*
 * Functions for computing the transformed distance of an edge (Na,Nb).
 */

/*
 * The C_EXPLICIT function returns the distance by looking it up in a table.
 */

int C_EXPLICIT(Node * Na, Node * Nb)
{
    return Na->Id < Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id];
}

/*
 * The C_FUNCTION function is used when the distance is defined by a
 * function (e.g. the Euclidean distance function). In order to speed
 * up the computations the following algorithm used:
 *
 *  (1) If (Na,Nb) is an edge on the current tour, then its distance 
 *      is available in either the field PredCost or SucCost.
 *
 *  (2) If the edge (Na,Nb) is a candidate edge incident to Na, then
 *      its distance is available in the field Cost of the corresponding
 *      Candidate structure.
 *     
 *  (3) A hash table (CacheVal) is consulted to see if the distance has
 *      been stored. 
 *      
 *  (4) Otherwise the distance function is called and the distance computed
 *      is stored in the hash table.                  
 */

int C_FUNCTION(Node * Na, Node * Nb)
{
    Node *Nc;
    Candidate *Cand;
    int Index, i, j;

    if (PredSucCostAvailable) {
        if (Na->Suc == Nb)
            return Na->SucCost;
        if (Na->Pred == Nb)
            return Na->PredCost;
    }
    if ((Cand = Na->CandidateSet))
        for (; (Nc = Cand->To); Cand++)
            if (Nc == Nb)
                return Cand->Cost;
    if ((Cand = Nb->CandidateSet))
        for (; (Nc = Cand->To); Cand++)
            if (Nc == Na)
                return Cand->Cost;
    if ((Cand = Na->BackboneCandidateSet))
        for (; (Nc = Cand->To); Cand++)
            if (Nc == Nb)
                return Cand->Cost;
    if ((Cand = Nb->BackboneCandidateSet))
        for (; (Nc = Cand->To); Cand++)
            if (Nc == Na)
                return Cand->Cost;
    if (CacheSig == 0)
        return D(Na, Nb);
    i = Na->Id;
    j = Nb->Id;
    if (i > j) {
        int k = i;
        i = j;
        j = k;
    }
    Index = ((i << 8) + i + j) & CacheMask;
    if (CacheSig[Index] == i)
        return CacheVal[Index];
    CacheSig[Index] = i;
    return (CacheVal[Index] = D(Na, Nb));
}

int D_EXPLICIT(Node * Na, Node * Nb)
{
    return (Na->Id <
            Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id]) + Na->Pi + Nb->Pi;
}

int D_FUNCTION(Node * Na, Node * Nb)
{
    return (Fixed(Na, Nb) ? 0 : Distance(Na, Nb) * Precision) + Na->Pi +
        Nb->Pi;
}//################################### C.c end ###################################

//################################### CandidateReport.c begin ###################################
#include "LKH.h"

/*
 * The CandidateReport function prints the minimum, average and maximum
 * number of candidates associated with a node.
 */

void CandidateReport()
{
    int Min = INT_MAX, Max = 0, Fixed = 0, Count;
    GainType Sum = 0, Cost = 0;
    Node *N;
    Candidate *NN;

    N = FirstNode;
    do {
        Count = 0;
        if (N->CandidateSet)
            for (NN = N->CandidateSet; NN->To; NN++)
                Count++;
        if (Count > Max)
            Max = Count;
        if (Count < Min)
            Min = Count;
        Sum += Count;
        if (N->FixedTo1 && N->Id < N->FixedTo1->Id) {
            Fixed++;
            Cost += Distance(N, N->FixedTo1);
        }
        if (N->FixedTo2 && N->Id < N->FixedTo2->Id) {
            Fixed++;
            Cost += Distance(N, N->FixedTo2);
        }
    }
    while ((N = N->Suc) != FirstNode);
    printff("Cand.min = %d, Cand.avg = %0.1f, Cand.max = %d\n",
            Min, (double) Sum / Dimension, Max);
    if (Fixed > 0)
        printff("Edges.fixed = %d [Cost = " GainFormat "]\n", Fixed, Cost);
    if (MergeTourFiles >= 1) {
        Count = 0;
        N = FirstNode;
        do
            if (IsCommonEdge(N, N->MergeSuc[0]))
                Count++;
        while ((N = N->Suc) != FirstNode);
        printff("Edges.common = %d\n", Count);
    }
}
//################################### CandidateReport.c end ###################################

//################################### ChooseInitialTour.c begin ###################################
#include "LKH.h"

/*
 * The ChooseInitialTour function generates a pseudo-random initial tour.
 * The algorithm constructs a tour as follows.
 *
 * First, a random node N is chosen.
 *
 * Then, as long as no all nodes have been chosen, choose the next node to
 * follow N in the tour, NextN, and set N equal to NextN.
 *
 * NextN is chosen as follows:
 *
 *  (A) If possible, and Trial = 1, choose NextN such that
 *      (N,NextN) is an edge of a given initial tour.
 *  (B) Otherwise, if  possible, choose NextN such that (N,NextN) is a
 *      fixed edge, or is common to two or more tours to be merged.
 *  (C) Otherwise, if possible, choose NextN so that (N,NextN) is a
 *      candidate edge, the alpha-value of (N,NextN) is zero, and (N,NextN)
 *      belongs to the current best or next best tour.
 *  (D) Otherwise, if possible, choose NextN such that (N,NextN) is a
 *      candidate edge.
 *  (E) Otherwise, choose NextN at random among those nodes not already
 *      chosen.
 *
 *  When more than one node may be chosen, the node is chosen at random
 *  among the alternatives (a one-way list of nodes).
 *
 *  The sequence of chosen nodes constitutes the initial tour.
 */

void ChooseInitialTour()
{
    Node *N, *NextN, *FirstAlternative, *Last;
    Candidate *NN;
    int Alternatives, Count, i;

    if (KickType > 0 && Kicks > 0 && Trial > 1) {
        for (Last = FirstNode; (N = Last->BestSuc) != FirstNode; Last = N)
            Follow(N, Last);
        for (i = 1; i <= Kicks; i++)
            KSwapKick(KickType);
        return;
    }

  Start:
    /* Mark all nodes as "not chosen" by setting their V field to zero */
    N = FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != FirstNode);
    Count = 0;

    /* Choose FirstNode without two incident fixed or common candidate edges */
    do {
        if (FixedOrCommonCandidates(N) < 2)
            break;
    }
    while ((N = N->Suc) != FirstNode);
    FirstNode = N;

    /* Move nodes with two incident fixed or common candidate edges in
       front of FirstNode */
    for (Last = FirstNode->Pred; N != Last; N = NextN) {
        NextN = N->Suc;
        if (FixedOrCommonCandidates(N) == 2)
            Follow(N, Last);
    }

    /* Mark FirstNode as chosen */
    FirstNode->V = 1;
    N = FirstNode;

    /* Loop as long as not all nodes have been chosen */
    while (N->Suc != FirstNode) {
        FirstAlternative = 0;
        Alternatives = 0;
        Count++;

        /* Case A */
        if (FirstNode->InitialSuc && Trial == 1 &&
            Count <= InitialTourFraction * Dimension) {
            for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                 if (!NextN->V && NextN == N->InitialSuc) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        /* Case B */
        if (Alternatives == 0) {
            for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                if (!NextN->V && Fixed(N, NextN)) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        if (Alternatives == 0 && MergeTourFiles > 1) {
            for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                if (!NextN->V && IsCommonEdge(N, NextN)) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        if (Alternatives == 0 && Trial > 1 &&
            ProblemType != HCP && ProblemType != HPP) {
            /* Case C */
            for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                if (!NextN->V && FixedOrCommonCandidates(NextN) < 2 &&
                    NN->Alpha == 0 && (InBestTour(N, NextN) ||
                                       InNextBestTour(N, NextN))) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        if (Alternatives == 0) {
            /* Case D */
            for (NN = N->CandidateSet; NN && (NextN = NN->To); NN++) {
                if (!NextN->V && FixedOrCommonCandidates(NextN) < 2) {
                    Alternatives++;
                    NextN->Next = FirstAlternative;
                    FirstAlternative = NextN;
                }
            }
        }
        if (Alternatives == 0) {
            /* Case E (actually not really a random choice) */
            NextN = N->Suc;
            while ((FixedOrCommonCandidates(NextN) == 2 ||
                    Forbidden(N, NextN)) && NextN->Suc != FirstNode)
                NextN = NextN->Suc;
            if (FixedOrCommonCandidates(NextN) == 2 || Forbidden(N, NextN)) {
                FirstNode = N;
                goto Start;
            }
        } else {
            NextN = FirstAlternative;
            if (Alternatives > 1) {
                /* Select NextN at random among the alternatives */
                i = Random() % Alternatives;
                while (i--)
                    NextN = NextN->Next;
            }
        }
        /* Include NextN as the successor of N */
        Follow(NextN, N);
        N = NextN;
        N->V = 1;
    }
    if (Forbidden(N, N->Suc)) {
        FirstNode = N;
        goto Start;
    }
    if (MaxTrials == 0) {
        GainType Cost = 0;
        N = FirstNode;
        do
            Cost += C(N, N->Suc) - N->Pi - N->Suc->Pi;
        while ((N = N->Suc) != FirstNode);
        Cost /= Precision;
        if (Cost < BetterCost) {
            BetterCost = Cost;
            RecordBetterTour();
        }
    }
}

//################################### ChooseInitialTour.c end ###################################

//################################### Connect.c begin ###################################
#include "LKH.h"

/*
 * Let T be a minimum spanning tree on the graph, and let N1 be a node of
 * degree one in T. The Connect function determines a shortest edge emanating
 * from N1, but not in T. At return, the Next field of N1 points to the end 
 * node of the edge, and its NextCost field contains the cost of the edge. 
 * However, the search for the shortest edge is stopped if an edge shorter 
 * than a specified threshold (Max) is found.
*/

void Connect(Node * N1, int Max, int Sparse)
{
    Node *N;
    Candidate *NN1;
    int d;

    N1->Next = 0;
    N1->NextCost = INT_MAX;
    if (!Sparse || N1->CandidateSet == 0 ||
        N1->CandidateSet[0].To == 0 || N1->CandidateSet[1].To == 0) {
        /* Find the requested edge in a dense graph */
        N = FirstNode;
        do {
            if (N == N1 || N == N1->Dad || N1 == N->Dad)
                continue;
            if (FixedOrCommon(N1, N)) {
                N1->NextCost = D(N1, N);
                N1->Next = N;
                return;
            }
            if (!N1->FixedTo2 && !N->FixedTo2 &&
                !Forbidden(N1, N) &&
                (!c || c(N1, N) < N1->NextCost) &&
                (d = D(N1, N)) < N1->NextCost) {
                N1->NextCost = d;
                if (d <= Max)
                    return;
                N1->Next = N;
            }
        }
        while ((N = N->Suc) != FirstNode);
    } else {
        /* Find the requested edge in a sparse graph */
        for (NN1 = N1->CandidateSet; (N = NN1->To); NN1++) {
            if (N == N1->Dad || N1 == N->Dad)
                continue;
            if (FixedOrCommon(N1, N)) {
                N1->NextCost = NN1->Cost + N1->Pi + N->Pi;
                N1->Next = N;
                return;
            }
            if (!N1->FixedTo2 && !N->FixedTo2 &&
                !Forbidden(N1, N) &&
                (d = NN1->Cost + N1->Pi + N->Pi) < N1->NextCost) {
                N1->NextCost = d;
                if (d <= Max)
                    return;
                N1->Next = N;
            }
        }
    }
}
//################################### Connect.c end ###################################

//################################### Create_POPMUSIC_CandidateSet.c begin ###################################
#include "LKH.h"

#define maxNeighbors POPMUSIC_MaxNeighbors
#define trials POPMUSIC_Trials
#define NB_RES POPMUSIC_Solutions
#define SAMPLE_SIZE POPMUSIC_SampleSize

#define d(a, b) (((a) == (b) ? 0 : D(a, b)\
                 - (a)->Pi - (b)->Pi) / Precision)
#define less(a, b, dmin) (!c || (c(a, b)\
                          - (a)->Pi - (b)->Pi) / Precision < dmin)

static void build_path(int n, int *path, int nb_clust);
static void fast_POPMUSIC(int n, int *path, int R);
static void swap(int *a, int *b);
static int unif(int low, int high);
static void shuffle(int n, int *path);
static GainType length_path(int n, int *path);
static void path_threeOpt(int N, int **D, int *best_sol,
                          GainType * best_cost);

static Node **node;
static Node **node_path;

static void optimize_path(int n, int *path);

/*
 * The Create_POPMUSIC_CandidateSet function creates for each node 
 * a candidate set of K edges using the POPMUSIC algorithm.
 *
 * The function is called from the CreateCandidateSet function.
 *
 * Programmed by Keld Helsgaun and Eric Taillard, 2018.
 *
 * References:
 *
 *     É. D. Taillard and K. Helsgaun,
 *     POPMUSIC for the Travelling Salesman Problem.
 *     European Journal of Operational Research, 272(2):420-429 (2019).
 *
 *     K. Helsgaun,
 *     Using POPMUSIC for Candidate Set Generation in the
 *     Lin-Kernighan-Helsgaun TSP Solver.
 *     Technical Report, Computer Science, Roskilde University, 2018.
 */

void Create_POPMUSIC_CandidateSet(int K)
{
    int n, i, no_res, setInitialSuc, d, deleted;
    int *solution;
    GainType cost, costSum = 0;
    GainType costMin = PLUS_INFINITY, costMax = MINUS_INFINITY;
    Node *N;
    double startTime, entryTime;
    int InitialTourAlgorithmSaved = InitialTourAlgorithm;

    entryTime = GetTime();
    if (TraceLevel >= 2)
        printff("Creating POPMUSIC candidate set ...\n");
    AddTourCandidates();
    if (MaxCandidates == 0) {
        N = FirstNode;
        do {
            if (!N->CandidateSet)
                eprintf("MAX_CANDIDATES = 0: No candidates");
        } while ((N = N->Suc) != FirstNode);
        if (TraceLevel >= 2)
            printff("done\n");
        return;
    }

    /* Create a tour containing all fixed or common edges */
    InitialTourAlgorithm = WALK;
    ChooseInitialTour();
    InitialTourAlgorithm = InitialTourAlgorithmSaved;
    /* N->V == 1 iff N is going to be deleted */
    N = FirstNode;
    do {
        N->V = FixedOrCommon(N, N->Pred);
    } while ((N = N->Suc) != FirstNode);

    n = Dimension;
    solution = (int *) malloc((n + 1) * sizeof(int));
    node = (Node **) malloc((n + 1) * sizeof(Node *));
    node_path = (Node **) malloc((n + 1) * sizeof(Node *));

    for (no_res = 1; no_res <= NB_RES; no_res++) {
        /* Create set of non-deleted nodes */
        n = deleted = 0;
        N = FirstNode;
        do {
            if (!N->V) {
                node[solution[n] = n] = N;
                n++;
            } else
                deleted++;
        } while ((N = N->Suc) != FirstNode);
        shuffle(n, solution);
        solution[n] = solution[0];
        node[n] = node[solution[0]];
        startTime = GetTime();
        build_path(n, solution, SAMPLE_SIZE);
        if (deleted > 0) {
            /* Create a one-way list representing the built tour */
            for (i = 1; i <= n; i++)
                node[solution[i - 1]]->Next = node[solution[i]];
            /* Insert the deleted nodes in the one-way list */
            N = node[solution[0]];
            do {
                if (!N->V && N->Suc->V) {
                    Node *OldNext = N->Next;
                    do
                        N = N->Next = N->Suc;
                    while (N->Suc->V);
                    N->Next = OldNext;
                }
            } while ((N = N->Suc) != node[solution[0]]);
            /* Convert the one-way list to a solution array */
            N = node[solution[0]];
            n = 0;
            do {
                node[solution[n] = n] = N;
                n++;
            } while ((N = N->Next) != node[solution[0]]);
        }
        solution[n] = solution[0];
        node[n] = node[solution[0]];
        cost = length_path(n, solution);
        if (TraceLevel >= 2) {
            printff("%d: Initial cost:  GainType, ", no_res, cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff("Gap = %0.2f%%, ",
                        100.0 * (cost - Optimum) / Optimum);
            printff("Time: %0.2f sec.\n", GetTime() - startTime);
        }
        startTime = GetTime();
        fast_POPMUSIC(n, solution, SAMPLE_SIZE * SAMPLE_SIZE);
        solution[n] = solution[0];
        node[n] = node[solution[0]];
        cost = length_path(n, solution);
        if (TraceLevel >= 2) {
            printff("%d: Improved cost: GainType, ", no_res, cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff("Gap = %0.2f%%, ",
                        100.0 * (cost - Optimum) / Optimum);
            printff("Time: %0.2f sec.\n", GetTime() - startTime);
        }
        costSum += cost;
        if (cost > costMax)
            costMax = cost;
        setInitialSuc = 0;
        if (cost < costMin) {
            costMin = cost;
            setInitialSuc = POPMUSIC_InitialTour && !InitialTourFileName;
        }
        for (i = 0; i < n; i++) {
            Node *a = node[solution[i]];
            Node *b = node[solution[i + 1]];
            d = D(a, b);
            AddCandidate(a, b, d, 1);
            AddCandidate(b, a, d, 1);
            if (setInitialSuc)
                a->InitialSuc = b;
        }
    }
    if (TraceLevel >= 2) {
        printff
            ("Cost.min = " GainFormat ", Cost.avg = %0.2f, Cost.max = "
             GainFormat "\n", costMin, (double) costSum / NB_RES, costMax);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff
                ("Gap.min = %0.2f%%, Gap.avg = %0.2f%%, Gap.max = %0.2f%%\n",
                 100.0 * (costMin - Optimum) / Optimum,
                 100.0 * ((double) costSum / NB_RES - Optimum) / Optimum,
                 100.0 * (costMax - Optimum) / Optimum);
    }
    free(solution);
    free(node);
    free(node_path);
    ResetCandidateSet();
    if (K > 0)
        TrimCandidateSet(K);
    AddTourCandidates();
    if (CandidateSetSymmetric)
        SymmetrizeCandidateSet();
    if (TraceLevel >= 2) {
        CandidateReport();
        printff("POPMUSIC Time = %0.2f sec.\n", GetTime() - entryTime);
        printff("done\n");
    }
}

/************************ Compute the length of a path ************************/
static GainType length_path(int n, int *path)
{
    Node *a, *b;
    GainType length = 0;
    int i;

    for (i = 1, a = node[path[0]]; i <= n; i++, a = b) {
        b = node[path[i]];
        length += d(a, b);
    }
    return length;
}

/*************** Optimize the path between path[0] and path[n] ****************/
static void optimize_path(int n, int *path)
{
    int i, j;
    int *order, **d;
    GainType length;
    Node *a, *b;

    order = (int *) malloc((n + 1) * sizeof(int));
    d = (int **) malloc((n + 1) * sizeof(int *));
    for (i = 0; i <= n; i++) {
        order[i] = i;
        d[i] = (int *) malloc((n + 1) * sizeof(int));
        node_path[i] = node[path[i]];
    }
    for (i = 0; i < n; i++) {
        a = node_path[i];
        for (j = i + 1; j <= n; j++) {
            b = node_path[j];
            d[i][j] = d[j][i] = 
                IsPossibleCandidate(a, b) ? d(a, b) : INT_MAX;
        }
    }
    d[0][n] = d[n][0] = 0;
    length = 0;
    for (i = 0; i < n; i++)
        length += d[i][i + 1];
    path_threeOpt(n, d, order, &length);
    for (i = 0; i <= n; i++)
        free(d[i]);
    free(d);
    for (i = 0; i <= n; i++)
        order[i] = path[order[i]];
    for (i = 0; i <= n; i++)
        path[i] = order[i];
    free(order);
}

/*********** Build recursively a path between path[0] and path[n] *************/
static void build_path(int n, int *path, int nb_clust)
{
    int i, j, d, dmin, closest, start, end;
    int *tmp_path, *sample, *assignment, *start_clust, *assigned;

    if (n <= 2)
        return;
    if (n <= nb_clust * nb_clust) {
        optimize_path(n, path);
        return;
    }

    /* Temporary path */
    tmp_path = (int *) malloc((n + 1) * sizeof(int));
    for (i = 0; i <= n; i++)
        tmp_path[i] = path[i];

    /* Let tmp_path[1] be the closest city to path[0] */
    dmin = d(node[tmp_path[1]], node[path[0]]);
    closest = 1;
    for (i = 2; i < n; i++) {
        if (less(node[tmp_path[i]], node[path[0]], dmin) &&
            (d = d(node[tmp_path[i]], node[path[0]])) < dmin) {
            dmin = d;
            closest = i;
        }
    }
    swap(tmp_path + 1, tmp_path + closest);

    /* Let tmp_path[2] be the closest city to path[n] */
    dmin = d(node[tmp_path[2]], node[path[n]]);
    closest = 2;
    for (i = 3; i < n; i++) {
        if (less(node[tmp_path[i]], node[path[n]], dmin) &&
            (d = d(node[tmp_path[i]], node[path[n]])) < dmin) {
            dmin = d;
            closest = i;
        }
    }
    swap(tmp_path + 2, tmp_path + closest);

    /* Choose a sample of nb_clust-2 random cities from tmp_path */
    for (i = 3; i <= nb_clust; i++)
        swap(tmp_path + i, tmp_path + unif(i, n - 1));
    swap(tmp_path + 2, tmp_path + nb_clust);
    sample = (int *) malloc((nb_clust + 2) * sizeof(int));
    for (i = 0; i <= nb_clust; i++)
        sample[i] = tmp_path[i];
    sample[nb_clust + 1] = path[n];
    optimize_path(nb_clust + 1, sample);

    /* Assign each city of path to the closest of the sample */
    assignment = (int *) malloc((n + 1) * sizeof(int));
    for (i = 1; i < n; i++) {
        dmin = INT_MAX;
        for (j = 1; j <= nb_clust; j++) {
            if (path[i] == sample[j]) {
                closest = j;
                break;
            }
            if (less(node[path[i]], node[sample[j]], dmin) &&
                (d = d(node[path[i]], node[sample[j]])) < dmin) {
                dmin = d;
                closest = j;
            }
        }
        assignment[i] = closest;
    }

    /* Build clusters: ith cluster has (start_clust[i+1]-start_clust[i])
       cities */
    start_clust = (int *) calloc(nb_clust + 1, sizeof(int));
    for (i = 1; i < n; i++)
        start_clust[assignment[i]]++;
    for (i = 1; i <= nb_clust; i++)
        start_clust[i] += start_clust[i - 1];

    /* Clusters are stored in tmp_path in the order given by sample */
    assigned = (int *) calloc(nb_clust + 1, sizeof(int));
    for (i = 1; i < n; i++)
        tmp_path[start_clust[assignment[i] - 1] +
                 assigned[assignment[i]]++] = path[i];

    /* Reorder original path */
    for (i = 1; i < n; i++)
        path[i] = tmp_path[i - 1];
    free(tmp_path);
    free(sample);
    free(assigned);
    free(assignment);
    for (i = 0; i < nb_clust; i++) {
        start = start_clust[i];
        end = start_clust[i + 1] + 1;
        /* Recursively optimize sub-path corresponding to each cluster */
        build_path(end - start, path + start, nb_clust);
    }
    free(start_clust);
}

static void reverse(int *path, int i, int j)
{
    while (i < j)
        swap(path + i++, path + j--);
}

static void circular_right_shift(int n, int *path, int positions)
{
    reverse(path, 0, positions - 1);
    reverse(path, positions, n - 1);
    reverse(path, 0, n - 1);
}

/* Fast POPMUSIC: Optimize independent subpaths of R cities */
static void fast_POPMUSIC(int n, int *path, int R)
{
    int scan, i;

    if (R > n)
        R = n;
    /* Optimize subpaths of R cities with R/2 overlap; 2 scans */
    for (scan = 1; scan <= 2; scan++) {
        if (scan == 2) {
            circular_right_shift(n, path, R / 2);
            path[n] = path[0];
            node[n] = node[path[0]];
        }
        for (i = 0; i < n / R; i++)
            optimize_path(R, path + R * i);
        if (n % R != 0) /* Optimize last portion of the path */
            optimize_path(R, path + n - R);
    }
}

/* Iterated 3-opt code */

static int n;
static int **dist;
static int *Create_POPMUSIC_CandidateSet_tour, *pos;
static int **neighbor;
static int *neighbors;
static int reversed;
static char *dontLook;
static GainType tourLength;

static void createNeighbors();
static void threeOpt();
static void doubleBridgeKick();
static int prev(int v);
static int next(int v);
static int PREV(int v);
static int NEXT(int v);

static void path_threeOpt(int N, int **D, int *best_sol,
                          GainType * best_cost)
{
    int i, j, trial;
    int *bestTour;
    GainType bestTourLength;

    n = N + 1;
    Create_POPMUSIC_CandidateSet_tour = (int *) malloc(n * sizeof(int));
    pos = (int *) malloc(n * sizeof(int));
    bestTour = (int *) malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        pos[bestTour[i] = Create_POPMUSIC_CandidateSet_tour[i] = best_sol[i]] = i;
    dist = D;
    createNeighbors();
    dontLook = (char *) calloc(n, sizeof(char));
    bestTourLength = tourLength = *best_cost;
    if (POPMUSIC_Trials == 0)
        trials = n;
    for (trial = 1; trial <= trials; trial++) {
        threeOpt();
        if (tourLength < bestTourLength) {
            for (i = 0; i < n; i++)
                bestTour[i] = Create_POPMUSIC_CandidateSet_tour[i];
            bestTourLength = tourLength;
        } else {
            for (i = 0; i < n; i++) {
                pos[Create_POPMUSIC_CandidateSet_tour[i] = bestTour[i]] = i;
                tourLength = bestTourLength;
            }
        }
        if (n <= 5 || trial == trials)
            break;
        doubleBridgeKick();
    }
    *best_cost = bestTourLength;
    for (i = 0; i < n; i++)
        pos[Create_POPMUSIC_CandidateSet_tour[i] = bestTour[i]] = i;
    reversed = next(0) == N;
    for (i = 0, j = 0; j < n; i = NEXT(i), j++)
        best_sol[j] = i;
    free(Create_POPMUSIC_CandidateSet_tour);
    free(pos);
    free(bestTour);
    free(neighbors);
    for (i = 0; i < n; i++)
        free(neighbor[i]);
    free(neighbor);
    free(dontLook);
}

static int unif(int low, int high)
{
    return low + Random() % (high - low + 1);
}

static void shuffle(int n, int *path)
{
    int i;
    for (i = 1; i < n; i++)
        swap(path + i, path + unif(0, i));
}

static void swap(int *a, int *b)
{
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

static int fixed(int a, int b)
{
    return (a == 0 && b == n - 1) || (a == n - 1 && b == 0) ||
        FixedOrCommon(node_path[a], node_path[b]);
}

static int prev(int v)
{
    return Create_POPMUSIC_CandidateSet_tour[pos[v] > 0 ? pos[v] - 1 : n - 1];
}

static int next(int v)
{
    return Create_POPMUSIC_CandidateSet_tour[pos[v] < n - 1 ? pos[v] + 1 : 0];
}

static int between(int v1, int v2, int v3)
{
    int a = pos[v1], b = pos[v2], c = pos[v3];
    return a <= c ? b >= a && b <= c : b <= c || b >= a;
}

static int PREV(int v)
{
    return reversed ? next(v) : prev(v);
}

static int NEXT(int v)
{
    return reversed ? prev(v) : next(v);
}

static int Create_POPMUSIC_CandidateSet_BETWEEN(int v1, int v2, int v3)
{
    return !reversed ? between(v1, v2, v3) : between(v3, v2, v1);
}

static void flip(int from, int to)
{
    int i, j, size, tmp;

    if (from == to)
        return;
    if (reversed) {
        tmp = from;
        from = to;
        to = tmp;
    }
    i = pos[from], j = pos[to];
    size = j - i;
    if (size < 0)
        size += n;
    if (size >= n / 2) {
        tmp = i;
        i = ++j < n ? j : 0;
        j = --tmp >= 0 ? tmp : n - 1;
    }
    while (i != j) {
        tmp = Create_POPMUSIC_CandidateSet_tour[i];
        pos[Create_POPMUSIC_CandidateSet_tour[i] = Create_POPMUSIC_CandidateSet_tour[j]] = i;
        pos[Create_POPMUSIC_CandidateSet_tour[j] = tmp] = j;
        if (++i == n)
            i = 0;
        if (i != j && --j < 0)
            j = n - 1;
    }
}

static void threeOpt()
{
    int improved = 1, a, b, c, d, e, f, xa, xc, xe, i, j;
    GainType g0, g1, g2, g3, gain;

    while (improved) {
        improved = 0;
        for (b = 0; b < n; b++) {
            if (dontLook[b])
                continue;
            dontLook[b] = 1;
            for (xa = 1; xa <= 2; xa++, reversed = !reversed) {
                a = PREV(b);
                if (fixed(a, b))
                    continue;
                g0 = dist[a][b];
                for (i = 0; i < neighbors[b]; i++) {
                    c = neighbor[b][i];
                    if (c == prev(b) || c == next(b))
                        continue;
                    g1 = g0 - dist[b][c];
                    if (g1 <= 0)
                        break;
                    for (xc = 1; xc <= 2; xc++) {
                        d = xc == 1 ? PREV(c) : NEXT(c);
                        if (d == a || fixed(c, d))
                            continue;
                        g2 = g1 + dist[c][d];
                        if (xc == 1) {
                            gain = g2 - dist[d][a];
                            if (gain > 0) {
                                flip(b, d);
                                tourLength -= gain;
                                dontLook[a] = dontLook[b] = 0;
                                dontLook[c] = dontLook[d] = 0;
                                improved = 1;
                                i = neighbors[b];
                                break;
                            }
                        }
                        for (j = 0; j < neighbors[d]; j++) {
                            e = neighbor[d][j];
                            if (e == prev(d) || e == next(d) ||
                                (xc == 2 && !Create_POPMUSIC_CandidateSet_BETWEEN(b, e, c)))
                                continue;
                            g3 = g2 - dist[d][e];
                            if (g3 <= 0)
                                break;
                            for (xe = 1; xe <= xc; xe++) {
                                if (xc == 1)
                                    f = Create_POPMUSIC_CandidateSet_BETWEEN(b, e,
                                                c) ? NEXT(e) : PREV(e);
                                else
                                    f = xe == 1 ? PREV(e) : NEXT(e);
                                if (f == a || fixed(e, f))
                                    continue;
                                gain = g3 + dist[e][f] - dist[f][a];
                                if (gain > 0) {
                                    if (xc == 1) {
                                        flip(b, d);
                                        if (f == PREV(e))
                                            flip(e, a);
                                        else
                                            flip(a, e);
                                    } else if (xe == 1) {
                                        flip(e, c);
                                        if (b == NEXT(a))
                                            flip(b, f);
                                        else
                                            flip(f, b);
                                    } else {
                                        flip(d, a);
                                        if (f == NEXT(e))
                                            flip(f, c);
                                        else
                                            flip(c, f);
                                        if (b == NEXT(d))
                                            flip(b, e);
                                        else
                                            flip(e, b);
                                    }
                                    tourLength -= gain;
                                    dontLook[a] = dontLook[b] = 0;
                                    dontLook[c] = dontLook[d] = 0;
                                    dontLook[e] = dontLook[f] = 0;
                                    improved = 1;
                                    goto Next_b;
                                }
                            }
                        }
                    }
                }
              Next_b:;
            }
        }
    }
}

static void createNeighbors()
{
    int i, j, k, d;

    neighbor = (int **) malloc(n * sizeof(int *));
    for (i = 0; i < n; i++)
        neighbor[i] = (int *) malloc((maxNeighbors + 1) * sizeof(int));
    neighbors = (int *) calloc(n, sizeof(int));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                continue;
            d = dist[i][j];
            k = neighbors[i] <
                maxNeighbors ? neighbors[i]++ : maxNeighbors;
            while (k > 0 && d < dist[i][neighbor[i][k - 1]]) {
                neighbor[i][k] = neighbor[i][k - 1];
                k--;
            }
            neighbor[i][k] = j;
        }
    }
}

static int legal(int r, int i, int t[4])
{
    while (--i >= 0)
        if (r == t[i])
            return 0;
    return !fixed(r, next(r));
}

static int select_t(int i, int t[4])
{
    int r, r0;
    r = r0 = unif(0, n - 1);
    while (!legal(r, i, t)) {
        if (++r == n)
            r = 0;
        if (r == r0)
            return -1;
    }
    return r;
}

static void doubleBridgeKick()
{
    int t[4], a, b, c, d, e, f, g, h, i;

    reversed = 0;
    for (i = 0; i <= 3; i++) {
        t[i] = select_t(i, t);
        if (t[i] < 0)
            return;
    }
    if (pos[t[0]] > pos[t[1]])
        swap(t, t + 1);
    if (pos[t[2]] > pos[t[3]])
        swap(t + 2, t + 3);
    if (pos[t[0]] > pos[t[2]])
        swap(t, t + 2);
    if (pos[t[1]] > pos[t[3]])
        swap(t + 1, t + 3);
    if (pos[t[1]] > pos[t[2]])
        swap(t + 1, t + 2);
    a = t[0];
    b = next(a);
    c = t[2];
    d = next(c);
    e = t[1];
    f = next(e);
    g = t[3];
    h = next(g);
    flip(b, c);
    if (f == next(e))
        flip(f, h);
    else
        flip(h, f);
    if (b == next(d))
        flip(c, d);
    else
        flip(d, c);
    dontLook[a] = dontLook[b] = dontLook[c] = dontLook[d] = 0;
    dontLook[e] = dontLook[f] = dontLook[g] = dontLook[h] = 0;
    tourLength -= dist[a][b] - dist[b][c] +
        dist[c][d] - dist[d][a] +
        dist[e][f] - dist[f][g] + dist[g][h] - dist[h][e];
}
//################################### Create_POPMUSIC_CandidateSet.c end ###################################

//################################### CreateCandidateSet.c begin ###################################
#include "LKH.h"

/*
 * The CreateCandidateSet function determines for each node its set of incident
 * candidate edges.
 *
 * The Ascent function is called to determine a lower bound on the optimal tour 
 * using subgradient optimization. But only if the penalties (the Pi-values) is
 * not available on file. In the latter case, the penalties is read from the 
 * file, and the lower bound is computed from a minimum 1-tree.      
 *
 * The function GenerateCandidates is called to compute the Alpha-values and to 
 * associate to each node a set of incident candidate edges.  
 *
 * The CreateCandidateSet function itself is called from LKHmain.
 */

void CreateCandidateSet()
{
    GainType Cost, MaxAlpha, A;
    Node *Na;
    int i;
    double EntryTime = GetTime();

    Norm = 9999;
    if (C == C_EXPLICIT) {
        Na = FirstNode;
        do {
            for (i = 1; i < Na->Id; i++)
                Na->C[i] *= Precision;
        }
        while ((Na = Na->Suc) != FirstNode);
    }
   
    if (TraceLevel >= 2)
        printff("Creating candidates ...\n");

    
    /* No PiFile specified or available */
    Na = FirstNode;
    do
        Na->Pi = 0;
    while ((Na = Na->Suc) != FirstNode);
    Cost = Ascent();

    LowerBound = (double) Cost / Precision;
    if (TraceLevel >= 1) {
        printff("Lower bound = %0.1f", LowerBound);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.2f%%",
                    100.0 * (Optimum - LowerBound) / Optimum);
        if (!PiFile)
            printff(", Ascent time = %0.2f sec.",
                    fabs(GetTime() - EntryTime));
        printff("\n");
    }
    MaxAlpha = (GainType) fabs(Excess * Cost);
    if ((A = Optimum * Precision - Cost) > 0 && A < MaxAlpha)
        MaxAlpha = A;
    if (CandidateSetType == POPMUSIC ||
        MaxCandidates == 0)
        OrderCandidateSet(MaxCandidates, MaxAlpha, CandidateSetSymmetric);
    else
        GenerateCandidates(MaxCandidates, MaxAlpha, CandidateSetSymmetric);

  End_CreateCandidateSet:
    ResetCandidateSet();
    if (MaxTrials > 0) {
        Na = FirstNode;
        do {
            if (!Na->CandidateSet || !Na->CandidateSet[0].To) {
                if (MaxCandidates == 0)
                    eprintf
                        ("MAX_CANDIDATES = 0: Node %d has no candidates",
                         Na->Id);
                else
                    eprintf("Node %d has no candidates", Na->Id);
            }
        }
        while ((Na = Na->Suc) != FirstNode);

    }
    if (C == C_EXPLICIT) {
        Na = FirstNode;
        do
            for (i = 1; i < Na->Id; i++)
                Na->C[i] += Na->Pi + NodeSet[i].Pi;
        while ((Na = Na->Suc) != FirstNode);
    }
    if (TraceLevel >= 1) {
        CandidateReport();
        printff("Preprocessing time = %0.2f sec.\n",
                fabs(GetTime() - EntryTime));
    }
}
//################################### CreateCandidateSet.c end ###################################

//################################### Distance.c begin ###################################
#include "LKH.h"

/*
 * Functions for computing distances (see TSPLIB).
 *
 * The appropriate function is referenced by the function pointer Distance.
 */

int Distance_ATSP(Node * Na, Node * Nb)
{
    int n = DimensionSaved;
    if ((Na->Id <= n) == (Nb->Id <= n))
        return M;
    if (abs(Na->Id - Nb->Id) == n)
        return 0;
    return Na->Id <= n ? Na->C[Nb->Id - n] : Nb->C[Na->Id - n];
}

int Distance_EXPLICIT(Node * Na, Node * Nb)
{
    return Na->Id < Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id];
}
//################################### Distance.c end ###################################

//################################### eprintf.c begin ###################################
#include "LKH.h"
#include <stdarg.h>

/* 
 * The eprintf function prints an error message and exits.
 */

void eprintf(const char *fmt, ...)
{
    va_list args;

    if (LastLine && *LastLine)
        fprintf(stderr, "\n%s\n", LastLine);
    fprintf(stderr, "\n*** Error ***\n");
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}
//################################### eprintf.c end ###################################

//################################### ERXT.c begin ###################################
#include "LKH.h"

/* 
 * The ERXT function applies the Edge Recombination Crossover operator (ERX)
 * on the two individuals (tours) represented by the Suc and Next references,
 * resepectively.
 * 
 * ERX was originally described in 
 *
 *     D. Whitley, T. Starkweather, and D. Fuquay,
 *     Scheduling Problems and the Traveling Salesman:
 *     the Genetic Edge Recombination Operator. 
 *     Proc. Third Int. Conf. on Genetic Algorithms and Their Applications
 *     (1989)
 * 
 * ERXT implements the variant of ERX based on tabu-edges (Edge-T) described in
 *
 *     Chuan-Kang Ting,
 *     Improving Edge Recombination through Alternate Inheritance and
 *     Greedy Manner.
 *     Lecture Notes in Computer Science 3004, pp. 207-216 (2004)
 *
 * However, ERXT does not implement the greedy strategy used by Edge-T for 
 * choosing foreign edges.
 */

static Node *FirstFree;
static int Tabu;

static Node *SelectNext(Node * N);

void ERXT()
{
    Node *N, *Next;
    int i;

    Tabu = 0;
    N = FirstNode;
    do {
        N->OldSuc = N->Suc;
        N->OldSuc->OldPred = N;
        N->Next->Prev = N;
        N->Suc->Pred = N;
        N->V = 0;
    }
    while ((N = N->Suc) != FirstNode);
    if (Dimension == DimensionSaved)
        FirstNode = &NodeSet[1 + Random() % Dimension];
    else {
        for (i = Random() % Dimension; i > 0; i--)
            FirstNode = FirstNode->Suc;
        if (FirstNode->Id <= DimensionSaved)
            FirstNode += DimensionSaved;
    }
    N = FirstNode;
    N->V = 1;
    FirstFree = N->Suc;
    N->Pred->Suc = N->Suc;
    N->Suc->Pred = N->Pred;
    for (i = 1; i < Dimension; i++) {
        Next = SelectNext(N);
        if (Next == FirstFree)
            FirstFree = Next->Suc;
        Next->Pred->Suc = Next->Suc;
        Next->Suc->Pred = Next->Pred;
        Link(N, Next);
        N = Next;
        N->V = 1;
    }
    Link(N, FirstNode);
}

/*
 * The EdgeCount function computes the number of unused edges emanating
 * from a given node, N.
 */

static int EdgeCount(Node * N)
{
    int Count = 0;
    Node *Next;

    if (!N->OldPred->V)
        Count++;
    if (!N->OldSuc->V)
        Count++;
    Next = N->Prev;
    if (!Next->V && Next != N->OldPred && Next != N->OldSuc)
        Count++;
    Next = N->Next;
    if (!Next->V && Next != N->OldPred && Next != N->OldSuc)
        Count++;
    return Count;
}

#define ERXT_IsCommonEdge(Na, Nb)\
    (((Na)->OldPred == (Nb) || (Na)->OldSuc == (Nb)) &&\
     ((Na)->Prev == (Nb) || (Na)->Next == (Nb)))

/*
 * The SelectNext function select the next node to be added as a neighbor
 * to a given node, N.
 *
 * The function chooses the neighbor node with the highest priority (See Ting's
 * paper). If two or more possible neighbor nodes have the same priority,
 * then one of them is chosen randomly. If the node has no neighbors on the two
 * two tours, then the first node in the list of unused nodes is chosen.
 */

static Node *SelectNext(Node * N)
{
    Node *Next, *Alternative[4];
    int Alternatives = 0, Score, MaxScore = INT_MIN, i;

    for (i = 1; i <= 4; i++) {
        Next = i == 1 ? N->OldPred : i == 2 ? N->OldSuc :
            i == 3 ? N->Prev : N->Next;
        if (!Next->V &&
            (i <= 2 || (Next != N->OldPred && Next != N->OldSuc))) {
            if (Fixed(N, Next) || ERXT_IsCommonEdge(N, Next))
                Score = INT_MAX;
            else
                Score = -EdgeCount(Next) - (i <= 2 ? Tabu : -Tabu);
            if (Score >= MaxScore) {
                if (Score > MaxScore)
                    Alternatives = 0;
                Alternative[Alternatives++] = Next;
                MaxScore = Score;
            }
        }
    }
    if (Alternatives > 0) {
        Next = Alternative[Random() % Alternatives];
        if (Next == N->OldPred || Next == N->OldSuc)
            Tabu++;
        if (Next == N->Prev || Next == N->Next)
            Tabu--;
        return Next;
    }
    Next = FirstFree;
    while (Forbidden(N, Next)) {
        Next = Next->Suc;
        if (Next == FirstFree)
            break;
    }
    return Next;
}
//################################### ERXT.c end ###################################

//################################### Excludable.c begin ###################################
#include "LKH.h"

/* 
 * The Excludable function is used to test if an edge, (ta,tb), 
 * of the tour may be excluded (when making a move). An edge is
 * excludable if it is on the original tour and has not previously 
 * been excluded (and inserted again) in the current series of moves.
 * If the edge is excludable, the function returns 1; otherwise 0.
 *
 * The function is called from the BestMove function in order to
 * test if the last edge to be excluded in a non-gainful r-opt move 
 * is excludable.  
 */

int Excludable(Node * ta, Node * tb)
{
    if (ta == tb->OldPred)
        return !tb->OldPredExcluded;
    if (ta == tb->OldSuc)
        return !tb->OldSucExcluded;
    return 0;
}
//################################### Excludable.c end ###################################

//################################### Exclude.c begin ###################################
#include "LKH.h"

/* 
 * The Exclude function is used to register that an edge, (ta,tb), 
 * of the original tour has been excluded in a move. This is done by
 * setting the appropriate flag, OldPredExcluded or OldSucExcluded, 
 * for each of the two end nodes.
 */

void Exclude(Node * ta, Node * tb)
{
    if (ta == tb->Pred || ta == tb->Suc)
        return;
    if (ta == tb->OldPred)
        tb->OldPredExcluded = 1;
    else if (ta == tb->OldSuc)
        tb->OldSucExcluded = 1;
    if (tb == ta->OldPred)
        ta->OldPredExcluded = 1;
    else if (tb == ta->OldSuc)
        ta->OldSucExcluded = 1;
}
//################################### Exclude.c end ###################################

//################################### FindTour.c begin ###################################
#include "LKH.h"

/*
 * After the candidate set has been created the FindTour function is called
 * a predetermined number of times (Runs).
 *
 * FindTour performs a number of trials, where in each trial it attempts
 * to improve a chosen initial tour using the modified Lin-Kernighan edge
 * exchange heuristics.
 *
 * Each time a better tour is found, the tour is recorded, and the candidates
 * are reorderded by the AdjustCandidateSet function. Precedence is given to
 * edges that are common to two currently best tours. The candidate set is
 * extended with those tour edges that are not present in the current set.
 * The original candidate set is re-established at exit from FindTour.
 */

static void SwapCandidateSets();
static GainType OrdinalTourCost;

GainType FindTour()
{
    GainType Cost;
    Node *t;
    int i;
    double EntryTime = GetTime();

    t = FirstNode;
    do
        t->OldPred = t->OldSuc = t->NextBestSuc = t->BestSuc = 0;
    while ((t = t->Suc) != FirstNode);
    if (Run == 1 && Dimension == DimensionSaved) {
        OrdinalTourCost = 0;
        for (i = 1; i < Dimension; i++)
            OrdinalTourCost += C(&NodeSet[i], &NodeSet[i + 1])
                - NodeSet[i].Pi - NodeSet[i + 1].Pi;
        OrdinalTourCost += C(&NodeSet[Dimension], &NodeSet[1])
            - NodeSet[Dimension].Pi - NodeSet[1].Pi;
        OrdinalTourCost /= Precision;
    }
    BetterCost = PLUS_INFINITY;
    if (MaxTrials > 0)
        HashInitialize(HTable);
    else {
        Trial = 1;
        ChooseInitialTour();
    }

    double recordTime = GetTime();
    GainType recordCost = BetterCost;
    for (Trial = 1; Trial <= MaxTrials; Trial++)
    {
        // if(Trial != 1){ // 保证至少跑完一次 Trial
        //     if (GetTime() - EntryTime >= TimeLimit ||
        //         GetTime() - StartTime >= TotalTimeLimit) {
        //         if (TraceLevel >= 1)
        //             printff("*** Time limit exceeded ***\n");
        //         Trial--;
        //         break;
        //     }
        //     // 耗尽了分配给此子问题的时间
        //     if(SubproblemSize > 0 && GetTime() - SubProblemStartTime >= SubProblemTotalTimeLimit){
        //         if (TraceLevel >= 1)
        //             printff("*** Time limit for this subproblem has been exhausted ***\n");
        //         Trial--;
        //         break;
        //     }
        //     // 在一定时间跨度内统计改进幅度
        //     if(GetTime() - recordTime >= TimeSpan){
        //         if(recordCost - BetterCost < TimeSpan*ScheduleScoreInSecond){
        //             if (TraceLevel >= 1)
        //                 printff("*** The extent of improvement("GainFormat") is too small in %.1fs ***\n",(recordCost - BetterCost), TimeSpan);
        //             Trial--;
        //             break;
        //         }
        //         else{
        //             recordTime = GetTime();
        //             recordCost = BetterCost;
        //         }
        //     }
        // }
        
        // 总时间限制
        if (GetTime() - StartTime >= TotalTimeLimit) {
                if (TraceLevel >= 1)
                    printff("*** Time limit exceeded in FindTour ***\n");
                Trial--;
                break;
        }

        // 在一定时间跨度内统计改进幅度
        if(GetTime() - recordTime >= TimeSpan){
            if(recordCost - BetterCost < TimeSpan*ScheduleScoreInSecond){
                if (TraceLevel >= 1)
                    printff("*** The extent of improvement("GainFormat") is too small in %.1fs ***\n",(recordCost - BetterCost), TimeSpan);
                Trial--;
                break;
            }
            else{
                recordTime = GetTime();
                recordCost = BetterCost;
            }
        }

        /* Choose FirstNode at random */
        if (Dimension == DimensionSaved)
            FirstNode = &NodeSet[1 + Random() % Dimension];
        else
            for (i = Random() % Dimension; i > 0; i--)
                FirstNode = FirstNode->Suc;
        ChooseInitialTour();
        Cost = LinKernighan();
        // if (Trial == 1 || (GetTime() - EntryTime < TimeLimit &&
        //                    GetTime() - StartTime < TotalTimeLimit)) {
        if (FirstNode->BestSuc) {
            /* Merge tour with current best tour */
            t = FirstNode;
            while ((t = t->Next = t->BestSuc) != FirstNode);
            Cost = MergeWithTour();
        }
        if (Dimension == DimensionSaved && Cost >= OrdinalTourCost &&
            BetterCost > OrdinalTourCost) {
            /* Merge tour with ordinal tour */
            for (i = 1; i < Dimension; i++)
                NodeSet[i].Next = &NodeSet[i + 1];
            NodeSet[Dimension].Next = &NodeSet[1];
            Cost = MergeWithTour();
        }
        // }
        if (Cost < BetterCost) {
            if (TraceLevel >= 1) {
                printff("* %d: Cost = " GainFormat, Trial, Cost);
                if (Optimum != MINUS_INFINITY && Optimum != 0)
                    printff(", Gap = %0.4f%%",
                            100.0 * (Cost - Optimum) / Optimum);
                printff(", Time = %0.2f sec. %s\n",
                        fabs(GetTime() - EntryTime),
                        Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
            }
            BetterCost = Cost;
            RecordBetterTour();

            if (StopAtOptimum && BetterCost == Optimum)
                break;
            AdjustCandidateSet();
            HashInitialize(HTable);
            HashInsert(HTable, Hash, Cost);
        } else if (TraceLevel >= 2)
            printff("  %d: Cost = " GainFormat ", Time = %0.2f sec.\n",
                    Trial, Cost, fabs(GetTime() - EntryTime));
        /* Record backbones if wanted */
        if (Trial <= BackboneTrials && BackboneTrials < MaxTrials) {
            SwapCandidateSets();
            AdjustCandidateSet();
            if (Trial == BackboneTrials) {
                if (TraceLevel >= 1) {
                    printff("# %d: Backbone candidates ->\n", Trial);
                    CandidateReport();
                }
            } else
                SwapCandidateSets();
        }
    }
    if (BackboneTrials > 0 && BackboneTrials < MaxTrials) {
        if (Trial > BackboneTrials ||
            (Trial == BackboneTrials &&
             (!StopAtOptimum || BetterCost != Optimum)))
            SwapCandidateSets();
        t = FirstNode;
        do {
            free(t->BackboneCandidateSet);
            t->BackboneCandidateSet = 0;
        } while ((t = t->Suc) != FirstNode);
    }
    t = FirstNode;
    if (Norm == 0 || MaxTrials == 0 || !t->BestSuc) {
        do
            t = t->BestSuc = t->Suc;
        while (t != FirstNode);
    }
    Hash = 0;
    do {
        (t->Suc = t->BestSuc)->Pred = t;
        Hash ^= Rand[t->Id] * Rand[t->Suc->Id];
    } while ((t = t->BestSuc) != FirstNode);
    if (Trial > MaxTrials)
        Trial = MaxTrials;
    ResetCandidateSet();
    return BetterCost;
}

/*
 * The SwapCandidateSets function swaps the normal and backbone candidate sets.
 */

static void SwapCandidateSets()
{
    Node *t = FirstNode;
    do {
        Candidate *Temp = t->CandidateSet;
        t->CandidateSet = t->BackboneCandidateSet;
        t->BackboneCandidateSet = Temp;
    } while ((t = t->Suc) != FirstNode);
}
//################################### FindTour.c end ###################################

//################################### FixedOrCommonCandidates.c begin ###################################
#include "LKH.h"

/* 
 * The FixedOrCommonCandidates function returns the number of fixed or
 * common candidate edges emanating from a given node, N.
 */

int FixedOrCommonCandidates(Node * N)
{
    int Count = 0;

    Count = N->FixedTo2 ? 2 : N->FixedTo1 ? 1 : 0;
    if (MergeTourFiles >= 2) {
        if (!Fixed(N, N->MergeSuc[0]) &&
            N->Subproblem == N->MergeSuc[0]->Subproblem &&
            IsCommonEdge(N, N->MergeSuc[0]))
            Count++;
        if (!Fixed(N->MergePred, N) &&
            N->Subproblem == N->MergePred->Subproblem &&
            IsCommonEdge(N->MergePred, N))
            Count++;
    }
    if (Count > 2)
        eprintf("Node %d has more than two required candidate edges",
                N->Id);
    return Count;
}
//################################### FixedOrCommonCandidates.c end ###################################

//################################### Flip.c begin ###################################
#include "LKH.h"

/*
 * The Flip function performs a 2-opt move. Edges (t1,t2) and (t3,t4) 
 * are exchanged with edges (t2,t3) and (t4,t1). Node t4 is one of 
 * t3's two neighbors on the tour; which one is uniquely determined
 * by the orientation of (t1,t2).
 *
 * The function is only used if the doubly linked list representation 
 * is used for a tour; if the two-level tree representation is used, 
 * the function Flip_SL is used instead.
 *
 * The 2-opt move is made by swapping Pred and Suc of each node of the
 * two segments, and then reconnecting the segments by suitable
 * settings of Pred and Suc of t1, t2, t3 and t4. In addition,
 * Rank is updated for nodes in the reversed segment (Rank gives the
 * ordinal number of a node in the tour).
 *
 * Any of two segments defined by the 2-opt move may be reversed. The
 * segment with the smallest number of nodes is reversed in order to
 * speed up computations. The number of nodes in a segment is found 
 * from the Rank-values. 
 * 
 * The move is pushed onto a stack of 2-opt moves. The stack makes it
 * possible to undo moves (by the RestoreTour function).
 *
 * Finally, the hash value corresponding to the tour is updated. 
 */

void Flip(Node * t1, Node * t2, Node * t3)
{
    Node *s1, *s2, *t4;
    int R, Temp, Ct2t3, Ct4t1;

    assert(t1->Pred == t2 || t1->Suc == t2);
    if (t3 == t2->Pred || t3 == t2->Suc)
        return;
    t4 = t1->Suc == t2 ? t3->Pred : t3->Suc;
    if (t1->Suc != t2) {
        s1 = t1;
        t1 = t2;
        t2 = s1;
        s1 = t3;
        t3 = t4;
        t4 = s1;
    }
    /* Find the segment with the smallest number of nodes */
    if ((R = t2->Rank - t3->Rank) < 0)
        R += Dimension;
    if (2 * R > Dimension) {
        s1 = t3;
        t3 = t2;
        t2 = s1;
        s1 = t4;
        t4 = t1;
        t1 = s1;
    }
    Ct2t3 = C(t2, t3);
    Ct4t1 = C(t4, t1);
    /* Swap segment (t3 --> t1) */
    R = t1->Rank;
    t1->Suc = 0;
    s2 = t3;
    while ((s1 = s2)) {
        s2 = s1->Suc;
        s1->Suc = s1->Pred;
        s1->Pred = s2;
        s1->Rank = R--;
        Temp = s1->SucCost;
        s1->SucCost = s1->PredCost;
        s1->PredCost = Temp;
    }
    (t3->Suc = t2)->Pred = t3;
    (t4->Suc = t1)->Pred = t4;
    t3->SucCost = t2->PredCost = Ct2t3;
    t1->PredCost = t4->SucCost = Ct4t1;
    SwapStack[Swaps].t1 = t1;
    SwapStack[Swaps].t2 = t2;
    SwapStack[Swaps].t3 = t3;
    SwapStack[Swaps].t4 = t4;
    Swaps++;
    Hash ^= (Rand[t1->Id] * Rand[t2->Id]) ^
        (Rand[t3->Id] * Rand[t4->Id]) ^
        (Rand[t2->Id] * Rand[t3->Id]) ^ (Rand[t4->Id] * Rand[t1->Id]);
}
//################################### Flip.c end ###################################

//################################### Flip_SL.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Flip_SL function performs a 2-opt move. Edges (t1,t2) and (t3,t4) 
 * are exchanged with edges (t2,t3) and (t4,t1). Node t4 is one of 
 * t3's two neighbors on the tour; which one is uniquely determined
 * by the orientation of (t1,t2).
 *
 * The function is only used if the two-level tree representation is used 
 * for a tour; if the doubly linked list representation is used, the function 
 * Flip is used instead.
 *
 * The worst-case time cost of a 2-op move is O(n) when the doubly linked
 * list representation is used. A worst-case cost of O(sqrt(n)) per 2-opt 
 * move may be achieved using the two-level tree representation.
 *
 * The idea is to divide the tour into roughly sqrt(n) segments. Each segment
 * is maintained as a doubly linked list of nodes (using pointers labeled
 * Pred and Suc). The segments are connected in a doubly linked list (using
 * pointers labeled Pred and Suc). Each segment contains a number, Rank,
 * that represents its position in the list, two pointers First and Last that
 * references the first and last node of the segment, respectively, and a bit,
 * Reversed, that is used to indicate whether the segment should be traversed
 * in forward or backward direction. Just switching this bit reverses the 
 * orientation of a whole segment. 
 *
 * The implementation of Flip_SL closely follows the suggestions given in
 *
 *      M. L. Fredman, D. S. Johnson & L. A. McGeoch,
 *      Data Structures for Traveling Salesmen",
 *      J. Algorithms, 16, 432-479 (1995).
 *
 * When a 2-opt move has been made it is pushed onto a stack of 2-opt moves.
 * The stack makes it possible to undo moves (by the RestoreTour function).
 *
 * Finally, the hash value corresponding to the tour is updated.
 */

static void Flip_SL_SplitSegment(Node * t1, Node * t2);

#define SPLIT_CUTOFF 0.75

void Flip_SL(Node * t1, Node * t2, Node * t3)
{
    Node *t4, *a, *b, *c, *d;
    Segment *P1, *P2, *P3, *P4, *Q1, *Q2;
    Node *s1, *s2;
    int i, Temp;

    assert(t1->Pred == t2 || t1->Suc == t2);
    if (t3 == t2->Pred || t3 == t2->Suc)
        return;
    if (Groups == 1) {
        Flip(t1, t2, t3);
        return;
    }
    t4 = t2 == SUC(t1) ? PRED(t3) : SUC(t3);
    P1 = t1->Parent;
    P2 = t2->Parent;
    P3 = t3->Parent;
    P4 = t4->Parent;
    /* Split segments if needed */
    if (P1 != P3 && P2 != P4) {
        if (P1 == P2) {
            Flip_SL_SplitSegment(t1, t2);
            P1 = t1->Parent;
            P2 = t2->Parent;
        }
        if (P3 == P4 && P1 != P3 && P2 != P4) {
            Flip_SL_SplitSegment(t3, t4);
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
    } else if ((P1 == P3
             && abs(t3->Rank - t1->Rank) > SPLIT_CUTOFF * GroupSize)
            || (P2 == P4
                && abs(t4->Rank - t2->Rank) > SPLIT_CUTOFF * GroupSize)) {
        if (P1 == P2) {
            Flip_SL_SplitSegment(t1, t2);
            P1 = t1->Parent;
            P2 = t2->Parent;
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
        if (P3 == P4) {
            Flip_SL_SplitSegment(t3, t4);
            P1 = t1->Parent;
            P2 = t2->Parent;
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
    }
    /* Check if it is possible to flip locally within a segment */
    b = 0;
    if (P1 == P3) {
        /* Either the t1 --> t3 path or the t2 --> t4 path lies 
           within one segment */
        if (t1->Rank < t3->Rank) {
            if (P1 == P2 && P1 == P4 && t2->Rank > t1->Rank) {
                a = t1;
                b = t2;
                c = t3;
                d = t4;
            } else {
                a = t2;
                b = t1;
                c = t4;
                d = t3;
            }
        } else {
            if (P1 == P2 && P1 == P4 && t2->Rank < t1->Rank) {
                a = t3;
                b = t4;
                c = t1;
                d = t2;
            } else {
                a = t4;
                b = t3;
                c = t2;
                d = t1;
            }
        }
    } else if (P2 == P4) {
        /* The t2 --> t4 path lies within one segment */
        if (t4->Rank < t2->Rank) {
            a = t3;
            b = t4;
            c = t1;
            d = t2;
        } else {
            a = t1;
            b = t2;
            c = t3;
            d = t4;
        }
    }
    if (b) {
        int Cbc = C(b, c), Cda = C(d, a);
        /* Flip locally (b --> d) within a segment */
        i = d->Rank;
        d->Suc = 0;
        s2 = b;
        while ((s1 = s2)) {
            s2 = s1->Suc;
            s1->Suc = s1->Pred;
            s1->Pred = s2;
            s1->Rank = i--;
            Temp = s1->SucCost;
            s1->SucCost = s1->PredCost;
            s1->PredCost = Temp;
        }
        d->Pred = a;
        b->Suc = c;
        d->PredCost = Cda;
        b->SucCost = Cbc;
        if (a->Suc == b) {
            a->Suc = d;
            a->SucCost = d->PredCost;
        } else {
            a->Pred = d;
            a->PredCost = d->PredCost;
        }
        if (c->Pred == d) {
            c->Pred = b;
            c->PredCost = b->SucCost;
        } else {
            c->Suc = b;
            c->SucCost = b->SucCost;
        }
        if (b->Parent->First == b)
            b->Parent->First = d;
        else if (d->Parent->First == d)
            d->Parent->First = b;
        if (b->Parent->Last == b)
            b->Parent->Last = d;
        else if (d->Parent->Last == d)
            d->Parent->Last = b;
    } else {
        int Ct2t3, Ct4t1;
        /* Reverse a sequence of segments */
        if (P1->Suc != P2) {
            a = t1;
            t1 = t2;
            t2 = a;
            a = t3;
            t3 = t4;
            t4 = a;
            Q1 = P1;
            P1 = P2;
            P2 = Q1;
            Q1 = P3;
            P3 = P4;
            P4 = Q1;
        }
        /* Find the sequence with the smallest number of segments */
        if ((i = P2->Rank - P3->Rank) < 0)
            i += Groups;
        if (2 * i > Groups) {
            a = t3;
            t3 = t2;
            t2 = a;
            a = t1;
            t1 = t4;
            t4 = a;
            Q1 = P3;
            P3 = P2;
            P2 = Q1;
            Q1 = P1;
            P1 = P4;
            P4 = Q1;
        }
        Ct2t3 = C(t2, t3);
        Ct4t1 = C(t4, t1);
        /* Reverse the sequence of segments (P3 --> P1). 
           Mirrors the corresponding code in the Flip function */
        i = P1->Rank;
        P1->Suc = 0;
        Q2 = P3;
        while ((Q1 = Q2)) {
            Q2 = Q1->Suc;
            Q1->Suc = Q1->Pred;
            Q1->Pred = Q2;
            Q1->Rank = i--;
            Q1->Reversed ^= 1;
        }
        P3->Suc = P2;
        P2->Pred = P3;
        P1->Pred = P4;
        P4->Suc = P1;
        if (t3->Suc == t4) {
            t3->Suc = t2;
            t3->SucCost = Ct2t3;
        } else {
            t3->Pred = t2;
            t3->PredCost = Ct2t3;
        }
        if (t2->Suc == t1) {
            t2->Suc = t3;
            t2->SucCost = Ct2t3;
        } else {
            t2->Pred = t3;
            t2->PredCost = Ct2t3;
        }
        if (t1->Pred == t2) {
            t1->Pred = t4;
            t1->PredCost = Ct4t1;
        } else {
            t1->Suc = t4;
            t1->SucCost = Ct4t1;
        }
        if (t4->Pred == t3) {
            t4->Pred = t1;
            t4->PredCost = Ct4t1;
        } else {
            t4->Suc = t1;
            t4->SucCost = Ct4t1;
        }
    }
    SwapStack[Swaps].t1 = t1;
    SwapStack[Swaps].t2 = t2;
    SwapStack[Swaps].t3 = t3;
    SwapStack[Swaps].t4 = t4;
    Swaps++;
    Hash ^= (Rand[t1->Id] * Rand[t2->Id]) ^
        (Rand[t3->Id] * Rand[t4->Id]) ^
        (Rand[t2->Id] * Rand[t3->Id]) ^ (Rand[t4->Id] * Rand[t1->Id]);
}

/*
   The Flip_SL_SplitSegment function is called by the Flip_SL function to split 
   a segment. Calling Flip_SL_SplitSegment(t1,t2), where t1 and t2 are neighbors 
   in the same segment, causes the segment to be split between t1 and t2. 
   The smaller half is merged with its neighbouring segment, thus keeping 
   the number of segments fixed.

   The implementation of Flip_SL_SplitSegment closely follows the suggestions given in

        M. L. Fredman, D. S. Johnson & L. A. McGeoch,
        Data Structures for Traveling Salesmen",
        J. Algorithms, 16, 432-479 (1995).
 */

void Flip_SL_SplitSegment(Node * t1, Node * t2)
{
    Segment *P = t1->Parent, *Q;
    Node *t, *u;
    int i, Temp, Count;

    if (t2->Rank < t1->Rank) {
        t = t1;
        t1 = t2;
        t2 = t;
    }
    Count = t1->Rank - P->First->Rank + 1;
    if (2 * Count < P->Size) {
        /* The left part of P is merged with its neighbouring segment, Q */
        Q = P->Reversed ? P->Suc : P->Pred;
        t = P->First->Pred;
        i = t->Rank;
        if (t == Q->Last) {
            if (t == Q->First && t->Suc != P->First) {
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Q->Reversed ^= 1;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            for (t = P->First; t != t2; t = t->Suc) {
                t->Parent = Q;
                t->Rank = ++i;
            }
            Q->Last = t1;
        } else {
            for (t = P->First; t != t2; t = u) {
                t->Parent = Q;
                t->Rank = --i;
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            Q->First = t1;
        }
        P->First = t2;
    } else {
        /* The right part of P is merged with its neighbouring segment, Q */
        Q = P->Reversed ? P->Pred : P->Suc;
        t = P->Last->Suc;
        i = t->Rank;
        if (t == Q->First) {
            if (t == Q->Last && t->Pred != P->Last) {
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Q->Reversed ^= 1;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            for (t = P->Last; t != t1; t = t->Pred) {
                t->Parent = Q;
                t->Rank = --i;
            }
            Q->First = t2;
        } else {
            for (t = P->Last; t != t1; t = u) {
                t->Parent = Q;
                t->Rank = ++i;
                u = t->Pred;
                t->Pred = t->Suc;
                t->Suc = u;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            Q->Last = t2;
        }
        Count = P->Size - Count;
        P->Last = t1;
    }
    P->Size -= Count;
    Q->Size += Count;
}
//################################### Flip_SL.c end ###################################

//################################### Flip_SSL.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Flip_SSL function performs a 2-opt move. Edges (t1,t2) and (t3,t4) 
 * are exchanged with edges (t2,t3) and (t4,t1). Node t4 is one of 
 * t3's two neighbors on the tour; which one is uniquely determined
 * by the orientation of (t1,t2).
 *
 * The function is only used if the three-level tree representation is used 
 * for a tour; if the doubly linked list representation is used, the Flip_SL
 * function is used instead.
 *
 * An average cost of O(n^(1/3)) per 2-opt move may be achieved using the 
 * three-level tree representation.
 *
 * See also the documentation for the Flip_SL function.
 */

static void Flip_SSL_SplitSegment(Node * t1, Node * t2);
static void SplitSSegment(Segment * t1, Segment * t2);
static void FlipNodes(Node * a, Node * b, Node * c, Node * d);
static void FlipSegments(Segment * a, Segment * b,
                         Segment * c, Segment * d);
static void FlipSSegments(SSegment * a, SSegment * b,
                          SSegment * c, SSegment * d);

#define SPLIT_CUTOFF 0.75

void Flip_SSL(Node * t1, Node * t2, Node * t3)
{
    Node *t4, *a, *b, *c, *d;
    Segment *P1, *P2, *P3, *P4;

    assert(t1->Pred == t2 || t1->Suc == t2);
    if (t3 == t2->Pred || t3 == t2->Suc)
        return;
    if (Groups == 1) {
        Flip(t1, t2, t3);
        return;
    }
    t4 = t2 == SUC(t1) ? PRED(t3) : SUC(t3);
    P1 = t1->Parent;
    P2 = t2->Parent;
    P3 = t3->Parent;
    P4 = t4->Parent;
    /* Split segments if needed */
    if (P1 != P3 && P2 != P4) {
        if (P1 == P2) {
            Flip_SSL_SplitSegment(t1, t2);
            P1 = t1->Parent;
            P2 = t2->Parent;
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
        if (P3 == P4 && P1 != P3 && P2 != P4) {
            Flip_SSL_SplitSegment(t3, t4);
            P1 = t1->Parent;
            P2 = t2->Parent;
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
    } else if ((P1 == P3
                && abs(t3->Rank - t1->Rank) > SPLIT_CUTOFF * GroupSize)
               || (P2 == P4
                   && abs(t4->Rank - t2->Rank) >
                   SPLIT_CUTOFF * GroupSize)) {
        if (P1 == P2) {
            Flip_SSL_SplitSegment(t1, t2);
            P1 = t1->Parent;
            P2 = t2->Parent;
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
        if (P3 == P4) {
            Flip_SSL_SplitSegment(t3, t4);
            P1 = t1->Parent;
            P2 = t2->Parent;
            P3 = t3->Parent;
            P4 = t4->Parent;
        }
    }
    /* Check if it is possible to flip locally within a segment */
    b = 0;
    if (P1 == P3) {
        /* Either the t1 --> t3 path or the t2 --> t4 path lies 
           within one segment */
        if (t1->Rank < t3->Rank) {
            if (P1 == P2 && P1 == P4 && t2->Rank > t1->Rank) {
                a = t1;
                b = t2;
                c = t3;
                d = t4;
            } else {
                a = t2;
                b = t1;
                c = t4;
                d = t3;
            }
        } else {
            if (P1 == P2 && P1 == P4 && t2->Rank < t1->Rank) {
                a = t3;
                b = t4;
                c = t1;
                d = t2;
            } else {
                a = t4;
                b = t3;
                c = t2;
                d = t1;
            }
        }
    } else if (P2 == P4) {
        /* The t2 --> t4 path lies within one segment */
        if (t4->Rank < t2->Rank) {
            a = t3;
            b = t4;
            c = t1;
            d = t2;
        } else {
            a = t1;
            b = t2;
            c = t3;
            d = t4;
        }
    }
    if (b)
        /* Flip locally (b --> d) within a segment */
        FlipNodes(a, b, c, d);
    else {
        Segment *a, *b, *c, *d, *t1 = P1, *t2 = P2, *t3 = P3, *t4 = P4;
        SSegment *P1, *P2, *P3, *P4, *Q1;
        P1 = t1->Parent;
        P2 = t2->Parent;
        P3 = t3->Parent;
        P4 = t4->Parent;
        if (P1 != P3 && P2 != P4) {
            if (P1 == P2) {
                SplitSSegment(t1, t2);
                P1 = t1->Parent;
                P2 = t2->Parent;
                P3 = t3->Parent;
                P4 = t4->Parent;
            }
            if (P3 == P4 && P1 != P3 && P2 != P4) {
                SplitSSegment(t3, t4);
                P1 = t1->Parent;
                P2 = t2->Parent;
                P3 = t3->Parent;
                P4 = t4->Parent;
            }
        } else
            if ((P1 == P3
                 && abs(t3->Rank - t1->Rank) > SPLIT_CUTOFF * SGroupSize)
                || (P2 == P4
                    && abs(t4->Rank - t2->Rank) >
                    SPLIT_CUTOFF * SGroupSize)) {
            if (P1 == P2) {
                SplitSSegment(t1, t2);
                P1 = t1->Parent;
                P2 = t2->Parent;
                P3 = t3->Parent;
                P4 = t4->Parent;
            }
            if (P3 == P4) {
                SplitSSegment(t3, t4);
                P1 = t1->Parent;
                P2 = t2->Parent;
                P3 = t3->Parent;
                P4 = t4->Parent;
            }
        }
        b = 0;
        if (P1 == P3) {
            if (t1->Rank < t3->Rank) {
                if (P1 == P2 && P1 == P4 && t2->Rank > t1->Rank) {
                    a = t1;
                    b = t2;
                    c = t3;
                    d = t4;
                } else {
                    a = t2;
                    b = t1;
                    c = t4;
                    d = t3;
                }
            } else {
                if (P1 == P2 && P1 == P4 && t2->Rank < t1->Rank &&
                    t1->Rank - t2->Rank + 1 < Groups) {
                    a = t3;
                    b = t4;
                    c = t1;
                    d = t2;
                } else {
                    a = t4;
                    b = t3;
                    c = t2;
                    d = t1;
                }
            }
        } else if (P2 == P4) {
            if (t4->Rank < t2->Rank) {
                a = t3;
                b = t4;
                c = t1;
                d = t2;
            } else {
                a = t1;
                b = t2;
                c = t3;
                d = t4;
            }
        }
        if (b)
            /* Flip locally (b --> d) within a super segment */
            FlipSegments(a, b, c, d);
        else {
            int i;
            if (P1->Suc != P2) {
                a = t1;
                t1 = t2;
                t2 = a;
                a = t3;
                t3 = t4;
                t4 = a;
                Q1 = P1;
                P1 = P2;
                P2 = Q1;
                Q1 = P3;
                P3 = P4;
                P4 = Q1;
            }
            /* Find the sequence with the fewest segments */
            if ((i = P2->Rank - P3->Rank) <= 0)
                i += SGroups;
            if (2 * i > SGroups) {
                a = t3;
                t3 = t2;
                t2 = a;
                a = t1;
                t1 = t4;
                t4 = a;
                Q1 = P3;
                P3 = P2;
                P2 = Q1;
                Q1 = P1;
                P1 = P4;
                P4 = Q1;
            }
            /* Reverse the sequence of segments (P3 --> P1) */
            FlipSSegments(P4, P3, P2, P1);
            if (t3->Suc == t4)
                t3->Suc = t2;
            else
                t3->Pred = t2;
            if (t2->Suc == t1)
                t2->Suc = t3;
            else
                t2->Pred = t3;
            if (t1->Pred == t2)
                t1->Pred = t4;
            else
                t1->Suc = t4;
            if (t4->Pred == t3)
                t4->Pred = t1;
            else
                t4->Suc = t1;
        }
    }
    if (!b) {
        int Ct2t3 = C(t2, t3), Ct4t1 = C(t4, t1);
        if (t3->Suc == t4) {
            t3->Suc = t2;
            t3->SucCost = Ct2t3;
        } else {
            t3->Pred = t2;
            t3->PredCost = Ct2t3;
        }
        if (t2->Suc == t1) {
            t2->Suc = t3;
            t2->SucCost = Ct2t3;
        } else {
            t2->Pred = t3;
            t2->PredCost = Ct2t3;
        }
        if (t1->Pred == t2) {
            t1->Pred = t4;
            t1->PredCost = Ct4t1;
        } else {
            t1->Suc = t4;
            t1->SucCost = Ct4t1;
        }
        if (t4->Pred == t3) {
            t4->Pred = t1;
            t4->PredCost = Ct4t1;
        } else {
            t4->Suc = t1;
            t4->SucCost = Ct4t1;
        }
    }
    SwapStack[Swaps].t1 = t1;
    SwapStack[Swaps].t2 = t2;
    SwapStack[Swaps].t3 = t3;
    SwapStack[Swaps].t4 = t4;
    Swaps++;
    Hash ^= (Rand[t1->Id] * Rand[t2->Id]) ^
        (Rand[t3->Id] * Rand[t4->Id]) ^
        (Rand[t2->Id] * Rand[t3->Id]) ^ (Rand[t4->Id] * Rand[t1->Id]);
}

/*
   The Flip_SSL_SplitSegment function is called by the Flip_SSL function to split a 
   segment. Calling Flip_SSL_SplitSegment(t1,t2), where t1 and t2 are neighbors in 
   the same segment, causes the segment to be split between t1 and t2. 
   The smaller half is merged with its neighboring segment, thus keeping
   the number of segments fixed.
 */

static void Flip_SSL_SplitSegment(Node * t1, Node * t2)
{
    Segment *P = t1->Parent, *Q;
    Node *t, *u;
    int i, Count, Temp;

    if (t2->Rank < t1->Rank) {
        t = t1;
        t1 = t2;
        t2 = t;
    }
    Count = t1->Rank - P->First->Rank + 1;
    if (2 * Count < P->Size) {
        /* The left part of P is merged with its neighbouring segment, Q */
        Q = P->Reversed ? P->Suc : P->Pred;
        t = P->First->Pred;
        i = t->Rank;
        if (t == Q->Last) {
            if (t == Q->First && t->Suc != P->First) {
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Q->Reversed ^= 1;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            for (t = P->First; t != t2; t = t->Suc) {
                t->Parent = Q;
                t->Rank = ++i;
            }
            Q->Last = t1;
        } else {
            for (t = P->First; t != t2; t = u) {
                t->Parent = Q;
                t->Rank = --i;
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            Q->First = t1;
        }
        P->First = t2;
    } else {
        /* The right part of P is merged with its neighbouring segment, Q */
        Q = P->Reversed ? P->Pred : P->Suc;
        t = P->Last->Suc;
        i = t->Rank;
        if (t == Q->First) {
            if (t == Q->Last && t->Pred != P->Last) {
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Q->Reversed ^= 1;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            for (t = P->Last; t != t1; t = t->Pred) {
                t->Parent = Q;
                t->Rank = --i;
            }
            Q->First = t2;
        } else {
            for (t = P->Last; t != t1; t = u) {
                t->Parent = Q;
                t->Rank = ++i;
                u = t->Pred;
                t->Pred = t->Suc;
                t->Suc = u;
                Temp = t->SucCost;
                t->SucCost = t->PredCost;
                t->PredCost = Temp;
            }
            Q->Last = t2;
        }
        Count = P->Size - Count;
        P->Last = t1;
    }
    P->Size -= Count;
    Q->Size += Count;
}

/*
   The SplitSSegment function is called by the Flip_SSL function to split a
   super segment. Calling SplitSSegment(t1,t2), where t1 and t2 are 
   neighbors in the same super segment, causes the super segment to be split 
   between t1 and t2.
   The smaller half is merged with its neighboring super segment, thus 
   keeping the number of super segments fixed.
 */

void SplitSSegment(Segment * t1, Segment * t2)
{
    SSegment *P = t1->Parent, *Q;
    Segment *t, *u;
    int i, Count;

    if (t2->Rank < t1->Rank) {
        t = t1;
        t1 = t2;
        t2 = t;
    }
    Count = t1->Rank - P->First->Rank + 1;
    if (2 * Count < P->Size) {
        /* The left part of P is merged with its neighbouring segment, Q */
        Q = P->Reversed ? P->Suc : P->Pred;
        t = P->First->Pred;
        i = t->Rank;
        if (t == Q->Last) {
            if (t == Q->First && t->Suc != P->First) {
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Q->Reversed ^= 1;
                t->Reversed ^= 1;
            }
            for (t = P->First; t != t2; t = t->Suc) {
                t->Parent = Q;
                t->Rank = ++i;
            }
            Q->Last = t1;
        } else {
            for (t = P->First; t != t2; t = u) {
                t->Parent = Q;
                t->Rank = --i;
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                t->Reversed ^= 1;
            }
            Q->First = t1;
        }
        P->First = t2;
    } else {
        /* The right part of P is merged with its neighbouring segment, Q */
        Q = P->Reversed ? P->Pred : P->Suc;
        t = P->Last->Suc;
        i = t->Rank;
        if (t == Q->First) {
            if (t == Q->Last && t->Pred != P->Last) {
                u = t->Suc;
                t->Suc = t->Pred;
                t->Pred = u;
                Q->Reversed ^= 1;
                t->Reversed ^= 1;
            }
            for (t = P->Last; t != t1; t = t->Pred) {
                t->Parent = Q;
                t->Rank = --i;
            }
            Q->First = t2;
        } else {
            for (t = P->Last; t != t1; t = u) {
                t->Parent = Q;
                t->Rank = ++i;
                u = t->Pred;
                t->Pred = t->Suc;
                t->Suc = u;
                t->Reversed ^= 1;
            }
            Q->Last = t2;
        }
        Count = P->Size - Count;
        P->Last = t1;
    }
    P->Size -= Count;
    Q->Size += Count;
}

static void FlipNodes(Node * a, Node * b, Node * c, Node * d)
{
    Node *s1, *s2 = b;
    int i = d->Rank, Temp, Cbc = C(b, c), Cda = C(d, a);

    d->Suc = 0;
    while ((s1 = s2)) {
        s2 = s1->Suc;
        s1->Suc = s1->Pred;
        s1->Pred = s2;
        s1->Rank = i--;
        Temp = s1->SucCost;
        s1->SucCost = s1->PredCost;
        s1->PredCost = Temp;
    }
    d->Pred = a;
    b->Suc = c;
    d->PredCost = Cda;
    b->SucCost = Cbc;
    if (a->Suc == b) {
        a->Suc = d;
        a->SucCost = d->PredCost;
    } else {
        a->Pred = d;
        a->PredCost = d->PredCost;
    }
    if (c->Pred == d) {
        c->Pred = b;
        c->PredCost = b->SucCost;
    } else {
        c->Suc = b;
        c->SucCost = b->SucCost;
    }
    if (b->Parent->First == b)
        b->Parent->First = d;
    else if (d->Parent->First == d)
        d->Parent->First = b;
    if (b->Parent->Last == b)
        b->Parent->Last = d;
    else if (d->Parent->Last == d)
        d->Parent->Last = b;
}

static void FlipSegments(Segment * a, Segment * b,
                         Segment * c, Segment * d)
{
    int i = d->Rank;
    Segment *s1, *s2 = b;

    d->Suc = 0;
    while ((s1 = s2)) {
        s2 = s1->Suc;
        s1->Suc = s1->Pred;
        s1->Pred = s2;
        s1->Rank = i--;
        s1->Reversed ^= 1;
    }
    d->Pred = a;
    b->Suc = c;
    if (a->Suc == b)
        a->Suc = d;
    else
        a->Pred = d;
    if (c->Pred == d)
        c->Pred = b;
    else
        c->Suc = b;
    if (b->Parent->First == b)
        b->Parent->First = d;
    else if (d->Parent->First == d)
        d->Parent->First = b;
    if (b->Parent->Last == b)
        b->Parent->Last = d;
    else if (d->Parent->Last == d)
        d->Parent->Last = b;
}

static void FlipSSegments(SSegment * a, SSegment * b,
                          SSegment * c, SSegment * d)
{
    int i = d->Rank;
    SSegment *s1, *s2 = b;

    d->Suc = 0;
    while ((s1 = s2)) {
        s2 = s1->Suc;
        s1->Suc = s1->Pred;
        s1->Pred = s2;
        s1->Rank = i--;
        s1->Reversed ^= 1;
    }
    d->Pred = a;
    b->Suc = c;
    if (a->Suc == b)
        a->Suc = d;
    else
        a->Pred = d;
    if (c->Pred == d)
        c->Pred = b;
    else
        c->Suc = b;
}
//################################### Flip_SSL.c end ###################################

//################################### Forbidden.c begin ###################################
#include "LKH.h"

/* 
 * The Forbidden function is used to test if an edge, (ta,tb), 
 * is one of the forbidden edges (C(ta, tb) == M) in a solution of
 * an asymmetric traveling saleman problem. 
 * If the edge is forbidden, the function returns 1; otherwise 0.
 */

int Forbidden(const Node * ta, const Node * tb)
{
    return ProblemType == ATSP &&
        (ta->Id <= DimensionSaved) == (tb->Id <= DimensionSaved);
}
//################################### Forbidden.c end ###################################

//################################### FreeStructures.c begin ###################################
#include "LKH.h"
#include "Sequence.h"
#include "Genetic.h"

/*      
 * The FreeStructures function frees all allocated structures.
 */

#define Free(s) { free(s); s = 0; }

void FreeStructures()
{
    FreeCandidateSets();
    FreeSegments();
    if (NodeSet) {
        int i;
        for (i = 1; i <= Dimension; i++) {
            Node *N = &NodeSet[i];
            Free(N->MergeSuc);
            N->C = 0;
        }
        Free(NodeSet);
    }
    Free(CostMatrix);
    Free(BestTour);
    Free(BetterTour);
    Free(SwapStack);
    Free(HTable);
    Free(Rand);
    Free(CacheSig);
    Free(CacheVal);
    Free(Name);
    Free(Type);
    Free(EdgeWeightType);
    Free(EdgeWeightFormat);
    Free(EdgeDataFormat);
    Free(NodeCoordType);
    Free(DisplayDataType);
    Free(Heap);
    Free(t);
    Free(T);
    Free(tSaved);
    Free(p);
    Free(q);
    Free(incl);
    Free(cycle);
    Free(G);
}

/*      
   The FreeSegments function frees the segments.
 */

void FreeSegments()
{
    if (FirstSegment) {
        Segment *S = FirstSegment, *SPrev;
        do {
            SPrev = S->Pred;
            Free(S);
        }
        while ((S = SPrev) != FirstSegment);
        FirstSegment = 0;
    }
    if (FirstSSegment) {
        SSegment *SS = FirstSSegment, *SSPrev;
        do {
            SSPrev = SS->Pred;
            Free(SS);
        }
        while ((SS = SSPrev) != FirstSSegment);
        FirstSSegment = 0;
    }
}

/*      
 * The FreeCandidateSets function frees the candidate sets.
 */

void FreeCandidateSets()
{
    Node *N = FirstNode;
    if (!N)
        return;
    do {
        Free(N->CandidateSet);
        Free(N->BackboneCandidateSet);
    }
    while ((N = N->Suc) != FirstNode);
}
//################################### FreeStructures.c end ###################################

//################################### Gain23.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Gain23 function attempts to improve a tour by making non-sequential
 * moves.
 *
 * The set of non-sequential moves considered consists of:
 *
 *  (1) any nonfeasible 2-opt move (producing two cycles) followed by any
 *      2- or 3-opt move which produces a feasible tour (by joining the
 *      two cycles);
 *
 *  (2) any nonfeasible 3-opt move (producing two cycles) followed by any
 *      2-opt move which produces a feasible tour (by joining the two
 *      cycles).
 *
 * The first and second move may in some cases be 4-opt in (2) and (1),
 * respectively. In (1) this can happen when a possible 3-opt move may
 * be extended to a nonfeasible 4-opt move. In (2) it can happen in those
 * cases where a sequential 4-opt extends a possible 3-opt move and
 * produces two cycles.
 *
 * The first move must have a positive gain. The second move is determined
 * by the BridgeGain function.
 */

/*
 The algorithm splits the set of possible moves up into a number disjoint
 subsets (called "cases"). When s1, s2, ..., s6 has been chosen, Case6 is
 used to discriminate between 7 cases. When s1, s2, ..., s8 has been
 chosen, Case8 is used to discriminate between 11 cases.
 
 A detailed description of the different cases can be found after the code.
 */

GainType Gain23()
{
    static Node *s1 = 0;
    static short OldReversed = 0;
    Node *s2, *s3, *s4, *s5, *s6 = 0, *s7, *s8 = 0, *s1Stop;
    Candidate *Ns2, *Ns4, *Ns6;
    GainType G0, G1, G2, G3, G4, G5, G6, Gain, Gain6;
    int X2, X4, X6, X8, Case6 = 0, Case8 = 0;
    int Breadth2, Breadth4, Breadth6;
    
    if (!s1 || s1->Subproblem != FirstNode->Subproblem)
        s1 = FirstNode;
    s1Stop = s1;
    for (X2 = 1; X2 <= 2; X2++) {
        Reversed = X2 == 1 ? OldReversed : (OldReversed ^= 1);
        do {
            s2 = SUC(s1);
            if (FixedOrCommon(s1, s2))
                continue;
            G0 = C(s1, s2);
            Breadth2 = 0;
            /* Choose (s2,s3) as a candidate edge emanating from s2 */
            for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
                if (s3 == s2->Pred || s3 == s2->Suc)
                    continue;
                if (++Breadth2 > MaxBreadth)
                    break;
                G1 = G0 - Ns2->Cost;
                for (X4 = 1; X4 <= 2; X4++) {
                    s4 = X4 == 1 ? SUC(s3) : PRED(s3);
                    if (FixedOrCommon(s3, s4))
                        continue;
                    G2 = G1 + C(s3, s4);
                    /* Try any gainful nonfeasible 2-opt move
                       followed by a 2-, 3- or 4-opt move */
                    if (X4 == 1 && s4 != s1 && !Forbidden(s4, s1) &&
                        2 * SegmentSize(s2, s3) <= Dimension &&
                        (!c || G2 - c(s4, s1) > 0) &&
                        (G3 = G2 - C(s4, s1)) > 0 &&
                        (Gain = BridgeGain(s1, s2, s3, s4, 0, 0, 0, 0, 0,
                                           G3)) > 0)
                        return Gain;
                    if (X4 == 2 &&
                        !Forbidden(s4, s1) &&
                        (!c || G2 - c(s4, s1) > 0) &&
                        (Gain = G2 - C(s4, s1)) > 0) {
                        Swap1(s1, s2, s3);
                        return Gain;
                    }
                    if (G2 - s4->Cost <= 0)
                        continue;
                    Breadth4 = 0;
                    /* Try any gainful nonfeasible 3- or 4-opt move
                       followed by a 2-opt move */
                    /* Choose (s4,s5) as a candidate edge emanating from s4 */
                    for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
                        if (s5 == s4->Pred || s5 == s4->Suc ||
                            (G3 = G2 - Ns4->Cost) <= 0)
                            continue;
                        if (++Breadth4 > MaxBreadth)
                            break;
                        /* Choose s6 as one of s5's two neighbors on the tour */
                        for (X6 = 1; X6 <= 2; X6++) {
                            if (X4 == 2) {
                                if (X6 == 1) {
                                    Case6 = 1 + !BETWEEN(s2, s5, s4);
                                    s6 = Case6 == 1 ? SUC(s5) : PRED(s5);
                                } else {
                                    s6 = s6 ==
                                    s5->Pred ? s5->Suc : s5->Pred;
                                    if (s5 == s1 || s6 == s1)
                                        continue;
                                    Case6 += 2;
                                }
                            } else if (BETWEEN(s2, s5, s3)) {
                                Case6 = 4 + X6;
                                s6 = X6 == 1 ? SUC(s5) : PRED(s5);
                                if (s6 == s1)
                                    continue;
                            } else {
                                if (X6 == 2)
                                    break;
                                Case6 = 7;
                                s6 = PRED(s5);
                            }
                            if (FixedOrCommon(s5, s6))
                                continue;
                            G4 = G3 + C(s5, s6);
                            Gain6 = 0;
                            if (!Forbidden(s6, s1) &&
                                s6 != s1->Pred &&
                                s6 != s1->Suc &&
                                (!c || G4 - c(s6, s1) > 0) &&
                                (Gain6 = G4 - C(s6, s1)) > 0) {
                                if (Case6 <= 2 || Case6 == 5 || Case6 == 6) {
                                    Make3OptMove(s1, s2, s3, s4, s5, s6,
                                                 Case6);
                                    return Gain6;
                                }
                                if ((Gain =
                                     BridgeGain(s1, s2, s3, s4, s5, s6, 0,
                                                0, Case6, Gain6)) > 0)
                                    return Gain;
                            }
                            Breadth6 = 0;
                            /* Choose (s6,s7) as a candidate edge
                               emanating from s6 */
                            for (Ns6 = s6->CandidateSet; (s7 = Ns6->To);
                                 Ns6++) {
                                if (s7 == s6->Pred || s7 == s6->Suc
                                    || (s6 == s2 && s7 == s3) || (s6 == s3
                                                                  && s7 ==
                                                                  s2)
                                    || (G5 = G4 - Ns6->Cost) <= 0)
                                    continue;
                                if (++Breadth6 > MaxBreadth)
                                    break;
                                /* Choose s8 as one of s7's two neighbors
                                 on the tour */
                                for (X8 = 1; X8 <= 2; X8++) {
                                    if (X8 == 1) {
                                        Case8 = Case6;
                                        switch (Case6) {
                                            case 1:
                                                s8 = BETWEEN(s2, s7,
                                                             s5) ? SUC(s7) :
                                                PRED(s7);
                                                break;
                                            case 2:
                                                s8 = BETWEEN(s3, s7,
                                                             s6) ? SUC(s7) :
                                                PRED(s7);
                                                break;
                                            case 3:
                                                if (BETWEEN(s5, s7, s4))
                                                    s8 = SUC(s7);
                                                else {
                                                    s8 = BETWEEN(s3, s7,
                                                                 s1) ? PRED(s7)
                                                    : SUC(s7);
                                                    Case8 = 17;
                                                }
                                                break;
                                            case 4:
                                                if (BETWEEN(s2, s7, s5))
                                                    s8 = BETWEEN(s2, s7,
                                                                 s4) ? SUC(s7)
                                                    : PRED(s7);
                                                else {
                                                    s8 = PRED(s7);
                                                    Case8 = 18;
                                                }
                                                break;
                                            case 5:
                                                s8 = PRED(s7);
                                                break;
                                            case 6:
                                                s8 = BETWEEN(s2, s7,
                                                             s3) ? SUC(s7) :
                                                PRED(s7);
                                                break;
                                            case 7:
                                                if (BETWEEN(s2, s7, s3))
                                                    s8 = SUC(s7);
                                                else {
                                                    s8 = BETWEEN(s5, s7,
                                                                 s1) ? PRED(s7)
                                                    : SUC(s7);
                                                    Case8 = 19;
                                                }
                                        }
                                    } else {
                                        if (Case8 >= 17 ||
                                            (Case6 != 3 && Case6 != 4
                                             && Case6 != 7))
                                            break;
                                        s8 = s8 ==
                                        s7->Pred ? s7->Suc : s7->Pred;
                                        Case8 += 8;
                                    }
                                    if (s8 == s1 ||
                                        (s7 == s1 && s8 == s2) ||
                                        (s7 == s3 && s8 == s4) ||
                                        (s7 == s4 && s8 == s3))
                                        continue;
                                    if (FixedOrCommon(s7, s8)
                                        || Forbidden(s8, s1))
                                        continue;
                                    G6 = G5 + C(s7, s8);
                                    if ((!c || G6 - c(s8, s1) > 0) &&
                                        (Gain = G6 - C(s8, s1)) > 0) {
                                        if (Case8 <= 15) {
                                            Make4OptMove(s1, s2, s3, s4,
                                                         s5, s6, s7, s8,
                                                         Case8);
                                            return Gain;
                                        }
                                        if (Gain > Gain6 &&
                                            (Gain =
                                             BridgeGain(s1, s2, s3, s4, s5,
                                                        s6, s7, s8, Case6,
                                                        Gain)) > 0)
                                            return Gain;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        while ((s1 = s2) != s1Stop);
    }
    return 0;
}

/*
 Below is shown the use of the variables X4, Case6 and Case8 to
 discriminate between cases considered by the algorithm.
 
 The notation
 
 ab-
 
 is used for a subtour that starts with the edge (sa,sb).
 For example, the tour
 
 12-43-
 
 contains the edges (s1,s2) and (s4,s3), in that order. A (*) signifies
 an infeasible solution. BridgeGain is called if the accumulated gain
 is possitive.
 
 X4 = 1:
 12-34-
 Case6 = 5:
 12-56-34- (*)
 Case8 = 5:
 12-87-56-34-, 12-56-87-34-, 12-56-34-87-
 Case8 = 13:
 12-78-56-34-, 12-56-78-34-, 12-56-34-78-
 Case6 = 6:
 12-65-34- (*)
 Case8 = 6:
 12-87-65-34-, 12-65-78-34-, 12-65-34-87-
 Case8 = 14:
 12-87-65-34-, 12-65-87-34-, 12-65-34-78-
 Case6 = 7:
 12-34-65- (*)
 Case8 = 7:
 12-78-34-65-
 Case8 = 15:
 12-87-34-65-
 Case8 = 19:
 12-34-87-65- (*), 12-34-78-65- (*)
 X4 = 2:
 12-43-
 Case6 = 1:
 12-56-43-
 Case8 = 1:
 12-78-56-43-, 12-56-87-43-, 12-56-43-87-
 Case6 = 2:
 12-43-65-
 Case8 = 2:
 12-87-43-65-, 12-43-78-65-, 12-43-65-87-
 Case6 = 3:
 12-65-43- (*)
 Case8 = 3:
 12-65-78-43-
 Case8 = 11:
 12-65-87-43-
 Case8 = 17:
 12-78-65-43- (*), 12-65-43-87- (*)
 Case6 = 4:
 12-43-56- (*)
 Case8 = 4:
 12-78-43-56-, 12-43-87-56-
 Case8 = 12:
 12-87-43-56-, 12-43-78-56-
 Case8 = 18:
 12-87-43-56- (*), 12-43-87-56- (*)
 */
//################################### Gain23.c end ###################################

//################################### GenerateCandidates.c begin ###################################
#include "LKH.h"

/*
 * The GenerateCandidates function associates to each node a set of incident 
 * candidate edges. The candidate edges of each node are sorted in increasing
 * order of their Alpha-values.
 *
 * The parameter MaxCandidates specifies the maximum number of candidate edges 
 * allowed for each node, and MaxAlpha puts an upper limit on their 
 * Alpha-values.
 *
 * A non-zero value of Symmetric specifies that the candidate set is to be
 * complemented such that every candidate edge is associated with both its 
 * two end nodes (in this way MaxCandidates may be exceeded). 
 *
 * The candidate edges of each node is kept in an array (CandidatSet) of
 * structures. Each structure (Candidate) holds the following information:
 *
 *      Node *To    : pointer to the other end node of the edge
 *      int Cost    : the cost (length) of the edge
 *      int Alpha   : the alpha-value of the edge
 *
 * The algorithm for computing Alpha-values in time O(n^2) and space O(n) 
 * follows the description in
 *
 *      Keld Helsgaun,
 *      An Effective Implementation of the Lin-Kernighan Traveling 
 *      Salesman Heuristic,
 *      Report, RUC, 1998. 
 */

static int Max(const int a, const int b)
{
    return a > b ? a : b;
}

void GenerateCandidates(int MaxCandidates, GainType MaxAlpha,
                        int Symmetric)
{
    Node *From, *To;
    Candidate *NFrom, *NN;
    int a, d, Count;

    if (TraceLevel >= 2)
        printff("Generating candidates ... ");
    if (MaxAlpha < 0 || MaxAlpha > INT_MAX)
        MaxAlpha = INT_MAX;
    /* Initialize CandidateSet for each node */
    FreeCandidateSets();
    From = FirstNode;
    do
        From->Mark = 0;
    while ((From = From->Suc) != FirstNode);

    if (MaxCandidates > 0) {
        do {
            From->CandidateSet =
               (Candidate *) malloc((MaxCandidates + 1) *
                                    sizeof(Candidate));
            From->CandidateSet[0].To = 0;
        }
        while ((From = From->Suc) != FirstNode);
    } else {
        AddTourCandidates();
        do {
            if (!From->CandidateSet)
                eprintf("MAX_CANDIDATES = 0: No candidates");
        } while ((From = From->Suc) != FirstNode);
        return;
    }

    /* Loop for each node, From */
    do {
        NFrom = From->CandidateSet;
        if (From != FirstNode) {
            From->Beta = INT_MIN;
            for (To = From; To->Dad != 0; To = To->Dad) {
                To->Dad->Beta =
                    !FixedOrCommon(To, To->Dad) ?
                    Max(To->Beta, To->Cost) : To->Beta;
                To->Dad->Mark = From;
            }
        }
        Count = 0;
        /* Loop for each node, To */
        To = FirstNode;
        do {
            if (To == From)
                continue;
            d = c && !FixedOrCommon(From, To) ? c(From, To) : D(From, To);
            if (From == FirstNode)
                a = To == From->Dad ? 0 : d - From->NextCost;
            else if (To == FirstNode)
                a = From == To->Dad ? 0 : d - To->NextCost;
            else {
                if (To->Mark != From)
                    To->Beta =
                        !FixedOrCommon(To, To->Dad) ?
                        Max(To->Dad->Beta, To->Cost) : To->Dad->Beta;
                a = d - To->Beta;
            }
            if (FixedOrCommon(From, To))
                a = INT_MIN;
            else {
                if (From->FixedTo2 || To->FixedTo2 || Forbidden(From, To))
                    continue;
                if (InInputTour(From, To)) {
                    a = 0;
                    if (c)
                        d = D(From, To);
                } else if (c) {
                    if (a > MaxAlpha ||
                        (Count == MaxCandidates &&
                         (a > (NFrom - 1)->Alpha ||
                          (a == (NFrom - 1)->Alpha
                           && d >= (NFrom - 1)->Cost))))
                        continue;
                    if (To == From->Dad) {
                        d = From->Cost;
                        a = 0;
                    } else if (From == To->Dad) {
                        d = To->Cost;
                        a = 0;
                    } else {
                        a -= d;
                        a += (d = D(From, To));
                    }
                }
            }
            if (a <= MaxAlpha && IsPossibleCandidate(From, To)) {
                /* Insert new candidate edge in From->CandidateSet */
                NN = NFrom;
                while (--NN >= From->CandidateSet) {
                    if (a > NN->Alpha || (a == NN->Alpha && d >= NN->Cost))
                        break;
                    *(NN + 1) = *NN;
                }
                NN++;
                NN->To = To;
                NN->Cost = d;
                NN->Alpha = a;
                if (Count < MaxCandidates) {
                    Count++;
                    NFrom++;
                }
                NFrom->To = 0;
            }
        }
        while ((To = To->Suc) != FirstNode);
    }
    while ((From = From->Suc) != FirstNode);

    AddTourCandidates();
    if (Symmetric)
        SymmetrizeCandidateSet();
    if (TraceLevel >= 2)
        printff("done\n");
}
//################################### GenerateCandidates.c end ###################################

//################################### GetTime.c begin ###################################
#define HAVE_GETRUSAGE
/* Undefine if you don't have the getrusage function */
/* #undef HAVE_GETRUSAGE */

/*
 * The GetTime function is used to measure execution time.
 *
 * The function is called before and after the code to be 
 * measured. The difference between the second and the
 * first call gives the number of seconds spent in executing
 * the code.
 *
 * If the system call getrusage() is supported, the difference 
 * gives the user time used; otherwise, the accounted real time.
 */

#ifdef HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>

// double GetTime() // 获取CPU时间
// {
//     struct rusage ru;
//     getrusage(RUSAGE_SELF, &ru);
//     return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1000000.0;
// }

double GetTime() // 获取实际时间
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}
#else

#include <time.h>

double GetTime()
{
    return (double) clock() / CLOCKS_PER_SEC;
}

#endif
//################################### GetTime.c end ###################################

//################################### gpx.c begin ###################################
/*****************************************************************************\
*                Generalized Partition Crossover 2                            *
*                                                                             *	
* Reference:  R.Tinos, D. Whitley, and G. Ochoa (2017).                       *
* A new generalized partition crossover for the traveling salesman problem:   *
* tunneling between local optima.                                             *
*                                                                             *
* Programmed by R. Tinos. Adapted for LKH by K. Helsgaun.                     *
*                                                                             *
\*****************************************************************************/
#include "LKH.h"
#include "gpx.h"

// GPX2
GainType gpx(int *solution_blue, int *solution_red, int *offspring)
{
    int i, j, *d4_vertices, n_d4_vertices, *common_edges_blue,
        *common_edges_red;
    int *common_edges_p2_blue, *common_edges_p2_red, *label_list,
        *label_list_inv, n_new;
    int *solution_blue_p2, *solution_red_p2, *vector_comp;
    GainType fitness_offspring;

    // Step 1: Identifying the vertices with degree 4 and the common edges
    d4_vertices = alloc_vectori(n_cities);
    common_edges_blue = alloc_vectori(n_cities);
    common_edges_red = alloc_vectori(n_cities);
    n_d4_vertices =
        d4_vertices_id(solution_blue, solution_red, d4_vertices,
                       common_edges_blue, common_edges_red);

    // Step 2: Insert ghost nodes
    n_new = n_cities + n_d4_vertices;
    // size of the new solutions: n_cities + number of ghost nodes
    label_list = alloc_vectori(n_new);
    // label_list: label for each node (including the ghost nodes)
    label_list_inv = alloc_vectori(n_cities);
    // label_list_inv: inverse of label_list
    j = 0;      // counter for the vertices with degree 4 (ghost nodes)
    for (i = 0; i < n_cities; i++) {
        label_list_inv[i] = -1;
        if (d4_vertices[i] == 1) {
            label_list[n_cities + j] = i;
            label_list_inv[i] = n_cities + j;
            j++;
        }
        label_list[i] = i;
    }
    // inserting the ghost nodes in solutions blue and red
    solution_blue_p2 = alloc_vectori(n_new);
    // solution blue with the ghost nodes
    solution_red_p2 = alloc_vectori(n_new);
    // solution red with the ghost nodes
    insert_ghost(solution_blue, solution_blue_p2, d4_vertices,
                 label_list_inv);
    insert_ghost(solution_red, solution_red_p2, d4_vertices,
                 label_list_inv);
    // identifying the common edges for the new solution
    common_edges_p2_blue = alloc_vectori(n_new);
    common_edges_p2_red = alloc_vectori(n_new);
    j = 0;
    for (i = 0; i < n_cities; i++) {
        common_edges_p2_blue[i] = common_edges_blue[i];
        common_edges_p2_red[i] = common_edges_red[i];
        if (d4_vertices[i] == 1) {
            common_edges_p2_blue[i] = 1;
            common_edges_p2_red[i] = 1;
            common_edges_p2_blue[n_cities + j] = common_edges_blue[i];
            common_edges_p2_red[n_cities + j] = common_edges_red[i];
            j++;
        }
    }

    // Step 3: creating the tour tables and finding the connected components
    vector_comp = alloc_vectori(n_new);
    // candidate component for each node (size n_new)
    // identify connected comp. using tour table
    tourTable(solution_blue_p2, solution_red_p2, solution_red, label_list,
              label_list_inv, vector_comp, n_new, common_edges_p2_blue,
              common_edges_p2_red);

    // Step 4: Creating the candidate components
    new_candidates(vector_comp, n_new); // object candidate recombination
    // component
    free(vector_comp);

    // Step 5: Finding the inputs and outputs of each candidate component
    findInputs(solution_blue_p2, solution_red_p2);

    // Step 6: testing the candidate components
    // Step 6.a: test components using simplified internal graphs
    for (i = 0; i < n_cand; i++)
        testComp(i);    // test component i

    // Step 6.b: test unfeasible components using simplified external graphs
    testUnfeasibleComp(solution_blue_p2);

    // Step 7.a: fusions of the candidate components that are neighbours
    //           (with more than to cutting points)
    fusion(solution_blue_p2, solution_red_p2);
    // if candidate i did not pass the test and has conditions, apply
    // fusion with the neighbour with more connections
    fusion(solution_blue_p2, solution_red_p2);
    // if candidate i did not pass the test and has conditions, apply
    // fusion with the neighbour with more connections
    fusion(solution_blue_p2, solution_red_p2);
    // if candidate i did not pass the test and has conditions, apply
    // fusion with the neighbour with more connections

    // Step 7.b: fusions of the candidate components in order to create
    //           partitions with two cutting points
    fusionB(solution_blue_p2, solution_red_p2);
    // if candidate i did not pass the test and has conditions, apply fusionB
    // to find fusions of partitions in order to have partitions with 2 cutting
    // points

    // Selecting the best between the blue and red path in each component
    fitness_offspring = off_gen(solution_blue_p2, solution_red_p2,
                                offspring, label_list);
    free_candidates();
    free(label_list);
    free(label_list_inv);
    free(d4_vertices);
    free(common_edges_blue);
    free(common_edges_p2_blue);
    free(common_edges_red);
    free(common_edges_p2_red);
    free(solution_blue_p2);
    free(solution_red_p2);
    return fitness_offspring;
}

// Compute the weights
int weight(int i, int j)
{
    Node *a = Map2Node[i], *b = Map2Node[j];
    return (C(a, b) - a->Pi - b->Pi) / Precision;
}

// Identifying the vertices with degree 4 and common edges
int d4_vertices_id(int *solution_blue, int *solution_red, int *d4_vertices,
                   int *common_edges_blue, int *common_edges_red)
{
    int i, aux, aux2, **M_aux, n_d4_vertices;

    //create a matrix (n_cities x 4) with all edges;
    M_aux = alloc_matrixi(n_cities, 4);

    for (i = 1; i < n_cities - 1; i++) {
        aux = solution_blue[i];
        M_aux[aux][0] = solution_blue[i + 1];
        M_aux[aux][1] = solution_blue[i - 1];
        aux = solution_red[i];
        M_aux[aux][2] = solution_red[i + 1];
        M_aux[aux][3] = solution_red[i - 1];
    }
    aux = solution_blue[0];
    M_aux[aux][0] = solution_blue[1];
    M_aux[aux][1] = solution_blue[n_cities - 1];
    aux = solution_red[0];
    M_aux[aux][2] = solution_red[1];
    M_aux[aux][3] = solution_red[n_cities - 1];
    aux = solution_blue[n_cities - 1];
    M_aux[aux][0] = solution_blue[0];
    M_aux[aux][1] = solution_blue[n_cities - 2];
    aux = solution_red[n_cities - 1];
    M_aux[aux][2] = solution_red[0];
    M_aux[aux][3] = solution_red[n_cities - 2];
    n_d4_vertices = 0; // number of degree 4 vertices

    for (i = 0; i < n_cities; i++) {
        d4_vertices[i] = 1;
        // d4_vertices: binary vector
        // (1: element is a degree 4 vertex; 0: otherwise);
        common_edges_blue[i] = 0;
        common_edges_red[i] = 0;
        aux = M_aux[i][0];
        aux2 = M_aux[i][2];
        if (aux == aux2 || aux == M_aux[i][3]) {
            d4_vertices[i] = 0;
            common_edges_blue[i] = 1;
            if (aux == aux2)
                common_edges_red[i] = 1;
        }
        aux = M_aux[i][1];
        if (aux == aux2 || aux == M_aux[i][3]) {
            d4_vertices[i] = 0;
            if (aux == aux2)
                common_edges_red[i] = 1;
        }
        if (d4_vertices[i] == 1)
            n_d4_vertices++;
    }

    dealloc_matrixi(M_aux, n_cities);
    return n_d4_vertices;
}

// Insert ghost nodes in the solution
void insert_ghost(int *solution, int *solution_p2, int *d4_vertices,
                  int *label_list_inv)
{
    int i, j, aux;

    j = 0;
    for (i = 0; i < n_cities; i++) {
        aux = solution[i];
        solution_p2[j] = aux;
        j++;
        if (d4_vertices[aux] == 1)
            solution_p2[j++] = label_list_inv[aux];
    }
}

// Finding the ghost pair (returns -1 if node has not a ghost pair)
int ghostPair(int *label_list, int *label_list_inv, int entry)
{
    return entry >
        n_cities - 1 ? label_list[entry] : label_list_inv[entry];
}

// Table code for the reverse solution
int tableCode(int ghost_a, int ghost_b, int ghost_c, int a, int b, int c,
              int common_a, int common_b, int ghost_flag)
{
    int ga, gb, gc;

    // vertices with degree 2
    if (common_a == 1 && common_b == 1)
        return -1;
    // vertices with degree 3 or 4
    if (ghost_a == -1)
        ga = 0;
    else
        ga = 1;
    if (ghost_b == -1)
        gb = 0;
    else
        gb = 1;
    if (ghost_c == -1)
        gc = 0;
    else
        gc = 1;
    if (ga == 0 && gb == 0 && gc == 0)
        return common_b == 1 ? a : c;
    if (ga == 0 && gb == 0 && gc == 1)
        return ghost_c;
    if (ga == 1 && gb == 0 && gc == 0)
        return a;
    if (ghost_flag == 0)
        return gc == 0 ? c : ghost_c;
    return a;
}

// correcting the number of entries 
// (removing common paths and assigned components)
// test if simplified graphs outside unfesible candidate component are equal
// Observation: this is equivalent of testing if all entries for a component
// are grouped after removing the feasible components (identified according
// to testComp) of the list of candidate entries
void simplifyPaths(int *solution_blue_p2, int n_new, int *vector_comp,
                   int *vector_cand, int *n_entries, int n_cand)
{
    int i, j, k, aux, *comp_seq, *inp_comp_seq;

    comp_seq = alloc_vectori(n_new);
    // sequence of components for all entries/exits in unfeasible components
    // in the order given by sol_blue
    inp_comp_seq = alloc_vectori(n_cand);
    // records the number of entries/exits in each component in comp_seq
    // creating comp_seq
    j = 0;      // j is the effective size of comp_seq
    k = solution_blue_p2[0];
    aux = vector_cand[k];
    if (vector_comp[k] == -1) {
        if (aux != vector_cand[solution_blue_p2[n_new - 1]])
            comp_seq[j++] = aux;
        if (aux != vector_cand[solution_blue_p2[1]])
            comp_seq[j++] = aux;
    }
    for (i = 1; i < n_new - 1; i++) {
        k = solution_blue_p2[i];
        aux = vector_cand[k];
        if (vector_comp[k] == -1) {
            if (aux != vector_cand[solution_blue_p2[i - 1]])
                comp_seq[j++] = aux;
            if (aux != vector_cand[solution_blue_p2[i + 1]])
                comp_seq[j++] = aux;
        }
    }
    k = solution_blue_p2[n_new - 1];
    aux = vector_cand[k];
    if (vector_comp[k] == -1) {
        if (aux != vector_cand[solution_blue_p2[n_new - 2]])
            comp_seq[j++] = aux;
        if (aux != vector_cand[solution_blue_p2[0]])
            comp_seq[j++] = aux;
    }
    for (i = 0; i < n_cand; i++)
        inp_comp_seq[i] = 0;

    // testing by checking the grouping of the components
    // (i.e., testing if the number of entries is 2)
    if (j > 0) {
        aux = comp_seq[0];
        if (aux != comp_seq[j - 1])
            inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        if (aux != comp_seq[1])
            inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        for (i = 1; i < j - 1; i++) {
            aux = comp_seq[i];
            if (aux != comp_seq[i - 1])
                inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
            if (aux != comp_seq[i + 1])
                inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        }
        aux = comp_seq[j - 1];
        if (aux != comp_seq[j - 2])
            inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        if (aux != comp_seq[0])
            inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        for (i = 0; i < n_cand; i++)
            if (n_entries[i] > 2 && inp_comp_seq[i] == 2)
                n_entries[i] = 2;
    }
    free(inp_comp_seq);
    free(comp_seq);
}

// Filling the first columns of the tour table
void tourTable_fill(int **Tour_table, int *d2_vertices,
                    int *solution_blue_p2, int *solution_red_p2,
                    int *solution_red, int *label_list,
                    int *label_list_inv, int *common_edges_blue_p2,
                    int *common_edges_red_p2, int n_new)
{
    int i, sol1, sol2, ghost_a, ghost_b, ghost_c, common_a, common_b,
        common_c, a, b, c;

    // Inserting in the table the blue and red tours (col. 0-1)
    for (i = 0; i < n_new - 1; i++) {
        // Inserting in the table the tour for the blue tour
        sol1 = solution_blue_p2[i];
        sol2 = solution_blue_p2[i + 1];
        if (common_edges_blue_p2[sol1] == 0) {
            Tour_table[sol1][0] = sol2;
            Tour_table[sol2][0] = sol1;
        } else {
            Tour_table[sol1][3] = sol2;
            Tour_table[sol2][3] = sol1;
        }
        // Inserting in the table the tour for the direct red tour
        sol1 = solution_red_p2[i];
        sol2 = solution_red_p2[i + 1];
        if (common_edges_red_p2[sol1] == 0) {
            Tour_table[sol1][1] = sol2;
            Tour_table[sol2][1] = sol1;
        }
    }
    sol1 = solution_blue_p2[n_new - 1];
    sol2 = solution_blue_p2[0];
    if (common_edges_blue_p2[sol1] == 0) {
        Tour_table[sol1][0] = sol2;
        Tour_table[sol2][0] = sol1;
    } else {
        Tour_table[sol1][3] = sol2;
        Tour_table[sol2][3] = sol1;
    }
    sol1 = solution_red_p2[n_new - 1];
    sol2 = solution_red_p2[0];
    if (common_edges_red_p2[sol1] == 0) {
        Tour_table[sol1][1] = sol2;
        Tour_table[sol2][1] = sol1;
    }
    // Inserting in the table the reverse red tours (col. 2)        
    a = solution_red[n_cities - 1];
    ghost_a = ghostPair(label_list, label_list_inv, a);
    if (ghost_a == -1)
        common_a = common_edges_red_p2[a];
    else
        common_a = common_edges_red_p2[ghost_a];
    b = solution_red[0];
    ghost_b = ghostPair(label_list, label_list_inv, b);
    if (ghost_b == -1)
        common_b = common_edges_red_p2[b];
    else
        common_b = common_edges_red_p2[ghost_b];
    c = solution_red[1];
    ghost_c = ghostPair(label_list, label_list_inv, c);
    if (ghost_c == -1)
        common_c = common_edges_red_p2[c];
    else
        common_c = common_edges_red_p2[ghost_c];
    Tour_table[b][2] = tableCode(ghost_a, ghost_b, ghost_c, a, b, c,
                                 common_a, common_b, 0);
    if (ghost_b != -1)
        Tour_table[ghost_b][2] = tableCode(ghost_a, ghost_b, ghost_c,
                                           a, b, c, common_a, common_b, 1);
    for (i = 1; i < n_cities - 1; i++) {
        a = b;
        ghost_a = ghost_b;
        common_a = common_b;
        b = c;
        ghost_b = ghost_c;
        common_b = common_c;
        c = solution_red[i + 1];
        ghost_c = ghostPair(label_list, label_list_inv, c);
        if (ghost_c == -1)
            common_c = common_edges_red_p2[c];
        else
            common_c = common_edges_red_p2[ghost_c];
        Tour_table[b][2] = tableCode(ghost_a, ghost_b, ghost_c,
                                     a, b, c, common_a, common_b, 0);
        if (ghost_b != -1)
            Tour_table[ghost_b][2] = tableCode(ghost_a, ghost_b, ghost_c,
                                               a, b, c, common_a, common_b,
                                               1);
    }
    a = b;
    ghost_a = ghost_b;
    common_a = common_b;
    b = c;
    ghost_b = ghost_c;
    common_b = common_c;
    c = solution_red[0];
    ghost_c = ghostPair(label_list, label_list_inv, c);
    if (ghost_c == -1)
        common_c = common_edges_red_p2[c];
    else
        common_c = common_edges_red_p2[ghost_c];
    Tour_table[b][2] = tableCode(ghost_a, ghost_b, ghost_c,
                                 a, b, c, common_a, common_b, 0);
    if (ghost_b != -1)
        Tour_table[ghost_b][2] = tableCode(ghost_a, ghost_b, ghost_c,
                                           a, b, c, common_a, common_b, 1);
}

// Identifying the vertices with degree 2
void d2_vertices_id(int *d2_vertices, int *solution_blue_p2,
                    int *common_edges_blue_p2, int n_new)
{
    int i;

    if (common_edges_blue_p2[solution_blue_p2[0]] == 1 &&
        common_edges_blue_p2[solution_blue_p2[n_new - 1]] == 1)
        d2_vertices[solution_blue_p2[0]] = 1;
    else
        d2_vertices[solution_blue_p2[0]] = 0;
    for (i = 1; i < n_new; i++) {
        if (common_edges_blue_p2[solution_blue_p2[i]] == 1 &&
            common_edges_blue_p2[solution_blue_p2[i - 1]] == 1)
            d2_vertices[solution_blue_p2[i]] = 1;
        else
            d2_vertices[solution_blue_p2[i]] = 0;
    }
}

// fixing the labels (in order to avoid gaps)
void labelsFix(int *vector_comp, int n_comp, int n_new)
{
    int i, j, *gap_labels, *new_label;

    gap_labels = alloc_vectori(n_comp);
    new_label = alloc_vectori(n_comp);
    for (i = 0; i < n_comp; i++) {
        gap_labels[i] = 1;
        new_label[i] = i;
    }
    for (i = 0; i < n_new; i++)
        gap_labels[vector_comp[i]] = 0;
    i = 0;
    j = n_comp - 1;
    while (j > i) {
        while (i < n_comp && gap_labels[i] == 0)
            i++;
        while (j >= 0 && gap_labels[j] == 1)
            j--;
        if (j > i) {
            new_label[j] = i;
            gap_labels[i] = 0;
            gap_labels[j] = 1;
        }
        i++;
        j--;
    }
    for (i = 0; i < n_new; i++)
        vector_comp[i] = new_label[vector_comp[i]];
    free(new_label);
    free(gap_labels);
}

// Finding the connected components using the tours table:
// one partition each time
void tourTable(int *solution_blue_p2, int *solution_red_p2,
               int *solution_red, int *label_list, int *label_list_inv,
               int *vector_comp, int n_new, int *common_edges_blue_p2,
               int *common_edges_red_p2)
{
    int i, k, cand_dir, cand_rev, n_comp, sol1, sol2, sol3, sol4, start,
        edge_tour, ghost_pair, n_rounds = 0, n_rounds_max = 1000;
    int min_size_dir, min_size_rev, red_chosen, cand_mcuts, cand_mcuts_dir,
        cand_mcuts_rev, min_size_dir_index, min_size_rev_index;
    int **Tour_table, *assigned_dir, *assigned_rev, *size_dir, *size_rev,
        *vector_comp_red;
    int *d2_vertices, *visited, *recently_assigned, *entries_flag_rev,
        *n_entries_dir, *n_entries_rev, *vector_cand_dir, *vector_cand_rev;

    // Memory allocation
    d2_vertices = alloc_vectori(n_new);
    visited = alloc_vectori(n_new);           // indicates the visited nodes
    vector_comp_red = alloc_vectori(n_new);   // indicates if ghost pair comes
                                              // from dir (0) or rev (1) red
    recently_assigned = alloc_vectori(n_new); // indicates the recently
                                              // assigned nodes
                                              // (for reversing ghost nodes)
    entries_flag_rev = alloc_vectori(n_new);  // auxiliary vector used for
                                              // checking direction of the
                                              // entries
    vector_cand_dir = alloc_vectori(n_new);   // auxiliary vector for
                                              // Tour_table (:,4)
    vector_cand_rev = alloc_vectori(n_new);   // auxiliary vector for
                                              // Tour_table (:,5)

    Tour_table = alloc_matrixi(n_new, 6);
    // Tours table
    // lines: vertices; 
    // columns: 
    //     0 - next single vertex in blue tour
    //     1 - next single vertex in direct red tour
    //     2 - next single vertex in reverse red tour
    //     3 - next common vertex 
    //     4 - candidate to connected component following the direct red tour
    //     5 - candidate to connected component following the reverse red tour
    // Obs.: all vertices has degree 3 or 2                                                                         
    // identifying the vertices with degree 2
    d2_vertices_id(d2_vertices, solution_blue_p2, common_edges_blue_p2,
                   n_new);
    // filling col. 0-3 of the tours table
    tourTable_fill(Tour_table, d2_vertices, solution_blue_p2,
                   solution_red_p2, solution_red, label_list,
                   label_list_inv, common_edges_blue_p2,
                   common_edges_red_p2, n_new);

    // remember that candidates with only one vertex should exist
    // (between common edges) connected components for vertices with degree 2
    // (each one has a label)
    n_comp = 0;
    for (i = 0; i < n_new; i++) {
        vector_comp_red[i] = -1;
        // -1 means that it was not assigned;
        // if assigned, can be 0 (dir. tour) or 1 (rev. tour)
        if (d2_vertices[i] == 1) {
            vector_comp[i] = n_comp;
            n_comp++;
        } else
            vector_comp[i] = -1; // indicates that vertex i was not
                                 // assigned yet                      
    }

    // finding the candidates to connected components (AB cycles) with any
    // number of cuts
    do {
        n_rounds++;
        // assigning the components 
        cand_mcuts_dir = 0;
        cand_mcuts_rev = 0;
        for (i = 0; i < n_new; i++) {
            if (vector_comp[i] == -1)
                visited[i] = 0;
            else
                visited[i] = 1;
            Tour_table[i][4] = -1;
            vector_cand_dir[i] = -1;
        }
        // folowing direct red tour
        // AB Cycles: direct red tour
        // all assigned become visited                          
        cand_dir = 0;
        for (i = 0; i < n_new; i++) {
            if (visited[i] == 0) {
                start = i;
                edge_tour = 0; // 0 for blue edge and 1 for red edge
                do {
                    Tour_table[i][4] = cand_dir;
                    vector_cand_dir[i] = cand_dir;
                    visited[i] = 1;
                    if (edge_tour == 0) {
                        i = Tour_table[i][0]; // get blue edge
                        edge_tour = 1;
                    } else {
                        i = Tour_table[i][1]; // get direct red edge
                        edge_tour = 0;
                    }
                } while (i != start);
                cand_dir++;
            }
        }
        // finding the number of entries and size       
        n_entries_dir = alloc_vectori(cand_dir);
        assigned_dir = alloc_vectori(cand_dir);
        size_dir = alloc_vectori(cand_dir);
        for (i = 0; i < cand_dir; i++) {
            n_entries_dir[i] = 0;
            size_dir[i] = 0;
        }
        for (i = 0; i < n_new; i++) {
            k = Tour_table[i][4];
            if (k != -1) {
                size_dir[k] = size_dir[k] + 1;
                if (k != Tour_table[Tour_table[i][3]][4]) {
                    n_entries_dir[k] = n_entries_dir[k] + 1;
                }
            }
        }
        // correcting the number of entries (removing common paths and
        // assigned components)
        simplifyPaths(solution_blue_p2, n_new, vector_comp,
                      vector_cand_dir, n_entries_dir, cand_dir);
        // following reverse red tour
        // AB Cycles: reverse red tour          
        // all assigned become visited
        cand_rev = 0;
        for (i = 0; i < n_new; i++) {
            if (vector_comp[i] == -1)
                visited[i] = 0;
            else
                visited[i] = 1;
            Tour_table[i][5] = -1;
            vector_cand_rev[i] = -1;
            entries_flag_rev[i] = 0;
            // 1 indicates that one of the entries for candidate cand_rev was
            // alredy assigned for direct red tour
            // obs.: the effective size is the number of candidates)
        }
        for (i = 0; i < n_new; i++) {
            if (visited[i] == 0) {
                start = i;
                edge_tour = 0; // 0 for blue edge and 1 for red edge
                do {
                    Tour_table[i][5] = cand_rev;
                    vector_cand_rev[i] = cand_rev;
                    ghost_pair = ghostPair(label_list, label_list_inv, i);
                    if (ghost_pair != -1) {
                        //check if i and ghost pair (if exists
                        if (vector_comp_red[i] == 0
                            || vector_comp_red[ghost_pair] == 0)
                            entries_flag_rev[cand_rev] = 1;
                        // 1 indicates that one of the entries for candidate 
                        // cand_rev was alredy assigned for direct red tour
                    }
                    visited[i] = 1;
                    if (edge_tour == 0) {
                        i = Tour_table[i][0]; // get blue edge
                        edge_tour = 1;
                    } else {
                        i = Tour_table[i][2]; // get reverse red edge
                        edge_tour = 0;
                    }
                } while (i != start);
                cand_rev++;
            }
        }
        // finding the number of entries and size
        n_entries_rev = alloc_vectori(cand_rev);
        assigned_rev = alloc_vectori(cand_rev);
        size_rev = alloc_vectori(cand_rev);
        for (i = 0; i < cand_rev; i++) {
            n_entries_rev[i] = 0;
            size_rev[i] = 0;
        }
        for (i = 0; i < n_new; i++) {
            k = Tour_table[i][5];
            if (k != -1) {
                size_rev[k] = size_rev[k] + 1;
                if (k != Tour_table[Tour_table[i][3]][5])
                    n_entries_rev[k] = n_entries_rev[k] + 1;
            }
        }
        // correcting the number of entries (removing common paths and
        // assigned components)
        simplifyPaths(solution_blue_p2, n_new, vector_comp,
                      vector_cand_rev, n_entries_rev, cand_rev);
        // Assigning the true candidates 
        // new labels for direct red tour
        min_size_dir = n_new; // minimum size for the candidates
        min_size_dir_index = -1;
        for (i = 0; i < cand_dir; i++) {
            cand_mcuts_dir++;
            assigned_dir[i] = n_comp; // new label
            n_comp++;
            if (size_dir[i] < min_size_dir ||
                (size_dir[i] == min_size_dir && n_entries_dir[i] == 2)) {
                min_size_dir = size_dir[i];
                min_size_dir_index = i;
            }
        }
        // new labels for reverse red tour
        min_size_rev = n_new; // minimum size for the candidates
        min_size_rev_index = -1;
        for (i = 0; i < cand_rev; i++) {
            if (entries_flag_rev[i] == 0) {
                cand_mcuts_rev++;
                assigned_rev[i] = n_comp; // new label
                n_comp++;
                if (size_rev[i] < min_size_rev ||
                    (size_rev[i] == min_size_rev
                     && n_entries_rev[i] == 2)) {
                    min_size_rev = size_rev[i];
                    min_size_rev_index = i;
                }
            } else
                assigned_rev[i] = -1;
        }
        cand_mcuts = cand_mcuts_dir + cand_mcuts_rev;
        if (cand_mcuts > 0 && n_rounds <= n_rounds_max) {
            // assigning components 
            // choose all components in one tour (only one) that has size
            // equal or smaller than the minimum size of the other component
            // use the number of entries when there is a tie
            if (min_size_rev < min_size_dir)
                red_chosen = 1; // 0 for direct and 1 for reverse
            else if (min_size_rev > min_size_dir)
                red_chosen = 0; // 0 for direct and 1 for reverse
            else {
                // Tie
                if (min_size_dir_index == -1)
                    red_chosen = 1;     // 0 for direct and 1 for reverse
                else {
                    if (min_size_rev_index == -1)
                        red_chosen = 0; // 0 for direct and 1 for reverse
                    else if (n_entries_dir[min_size_dir_index] == 2)
                        red_chosen = 0; // 0 for direct and 1 for reverse                                 
                    else if (n_entries_rev[min_size_rev_index] == 2)
                        red_chosen = 1; // 0 for direct and 1 for reverse                                 
                    else if (n_entries_rev[min_size_rev_index] <
                             n_entries_dir[min_size_dir_index])
                        red_chosen = 1; // 0 for direct and 1 for reverse 
                    else
                        red_chosen = 0; // 0 for direct and 1 for reverse
                }

            }
            for (i = 0; i < cand_rev; i++) {
                if (assigned_rev[i] != -1) {
                    if (red_chosen == 0 || size_rev[i] > min_size_dir) {
                        assigned_rev[i] = -1;
                        cand_mcuts_rev--;
                    }
                }
            }
            if (red_chosen == 1 && cand_mcuts_rev == 0) {
                red_chosen = 0;
                min_size_rev = min_size_dir;
            }
            for (i = 0; i < cand_dir; i++) {
                if (red_chosen == 1 || size_dir[i] > min_size_rev) {
                    assigned_dir[i] = -1;
                    cand_mcuts_dir--;
                }
            }
            cand_mcuts = cand_mcuts_dir + cand_mcuts_rev;
            if (cand_mcuts > 0) {
                for (i = 0; i < n_new; i++) {
                    if (vector_comp[i] == -1) {
                        if (red_chosen == 0
                            && assigned_dir[Tour_table[i][4]] != -1) {
                            vector_comp[i] =
                                assigned_dir[Tour_table[i][4]];
                            // assigning component
                            // recording dir. red tour                      
                            if (vector_comp_red[i] == -1) {
                                ghost_pair =
                                    ghostPair(label_list, label_list_inv,
                                              i);
                                if (ghost_pair != -1) {
                                    vector_comp_red[i] = 0;
                                    // zero means that it comes from the
                                    // direct red tour
                                    vector_comp_red[ghost_pair] = 0;
                                    // zero means that it comes from the
                                    // direct red tour
                                    // reversing the ghost nodes (changingi
                                    // direction in  Table) of the red tours
                                    // exchanging the reverse red edge for i
                                    // and ghost node
                                    sol1 = i;
                                    sol2 = ghost_pair;
                                    sol3 = Tour_table[sol1][2];
                                    sol4 = Tour_table[sol2][2];
                                    Tour_table[sol1][2] = sol4;
                                    Tour_table[sol4][2] = sol1;
                                    Tour_table[sol2][2] = sol3;
                                    Tour_table[sol3][2] = sol2;
                                }
                            }
                        } else if (red_chosen == 1
                                   && assigned_rev[Tour_table[i][5]] !=
                                   -1) {
                            vector_comp[i] =
                                assigned_rev[Tour_table[i][5]];
                            // assigning component                      
                            // recording dir. red tour                      
                            if (vector_comp_red[i] == -1) {
                                ghost_pair =
                                    ghostPair(label_list, label_list_inv,
                                              i);
                                if (ghost_pair != -1) {
                                    vector_comp_red[i] = 1;
                                    // zero means that it comes from the
                                    // direct red tour
                                    vector_comp_red[ghost_pair] = 1;
                                    // zero means that it comes from thei
                                    // direct red tour
                                    // reversing the ghost nodes (changingi
                                    // direction in Table) of the red tours
                                    // exchanging the reverse red edge for i
                                    // and ghost node
                                    sol1 = i;
                                    sol2 = ghost_pair;
                                    sol3 = Tour_table[sol1][1];
                                    sol4 = Tour_table[sol2][1];
                                    Tour_table[sol1][1] = sol4;
                                    Tour_table[sol4][1] = sol1;
                                    Tour_table[sol2][1] = sol3;
                                    Tour_table[sol3][1] = sol2;
                                }
                            }
                        }
                    }
                }
            }
        } else {
            // When the maximum number of rounds is reached, assign from
            // direct red tour
            for (i = 0; i < cand_dir; i++) {
                assigned_dir[i] = n_comp;
                n_comp++;
            }
            // assigning new labels 
            for (i = 0; i < n_new; i++) {
                if (vector_comp[i] == -1)
                    vector_comp[i] = assigned_dir[Tour_table[i][4]];
            }
        }
        free(n_entries_dir);
        free(assigned_dir);
        free(size_dir);
        free(n_entries_rev);
        free(assigned_rev);
        free(size_rev);
    } while (cand_mcuts > 0 && n_rounds <= n_rounds_max);

    labelsFix(vector_comp, n_comp, n_new); //fixing the labels

    // change ghost nodes in red tour
    for (i = 0; i < n_new; i++) {
        ghost_pair =
            ghostPair(label_list, label_list_inv, solution_red_p2[i]);
        if (vector_comp_red[solution_red_p2[i]] == 1 && ghost_pair > -1)
            solution_red_p2[i] = ghost_pair;
    }

    // Deallocating memory
    free(vector_comp_red);
    dealloc_matrixi(Tour_table, n_new);
    free(d2_vertices);
    free(visited);
    free(recently_assigned);
    free(entries_flag_rev);
    free(vector_cand_dir);
    free(vector_cand_rev);
}

/******************************************************************************\
*                                   Tour                                       *
\******************************************************************************/

static int *id;        // vector with the id(candidate component) of each node
static int *size;      // vector
static int *n_inputs;
static int *n_outputs;
static tour *blue, *red;
static unsigned int gpx_n; // size of the tours 
static int **M_neigh;  // neighbourhood matrix: first collumn indicates the
                       // number of neighbours, and the collumns 2 and 3
                       // indicate the index of the neighbours
static int **M_neigh2; // neighbourhood matrix: the collumns indicate the
                       // number of i conections to the neighbours indicated
                       // in collumns 2 and 3
int *test; // test of the candidates: 1 - true component; 0 - otherwise
static int isequal(Graph * G1, Graph * G2); // test if two graphs are equal

void new_candidates(int *vector_comp, int n_new)
{
    int i, j;

    gpx_n = n_new;
    n_cand = 0; // number of candidate components
    for (i = 0; i < gpx_n; i++)
        if (vector_comp[i] > n_cand)
            n_cand = vector_comp[i]; // remember that the first component has
                                     // label 0
    n_cand++;
    size = new_int(n_cand); // size of each component
    id = new_int(gpx_n);
    n_inputs = new_int(n_cand);
    n_outputs = new_int(n_cand);
    M_neigh = alloc_matrixi(n_cand, 3);
    M_neigh2 = alloc_matrixi(n_cand, 2);
    for (i = 0; i < n_cand; i++)
        size[i] = 0;
    for (i = 0; i < gpx_n; i++) {
        id[i] = vector_comp[i];
        j = id[i];
        size[j] = size[j] + 1;
    }
    test = new_int(n_cand);
    blue = new_tour(n_cand);
    red = new_tour(n_cand);
}

void free_candidates(void)
{
    int i;

    free(size);
    free(id);
    free(n_inputs);
    free(n_outputs);
    dealloc_matrixi(M_neigh, n_cand);
    dealloc_matrixi(M_neigh2, n_cand);
    for (i = 0; i < n_cand; i++) {
        free(blue[i].inputs);
        free(blue[i].outputs);
        free(red[i].inputs);
        free(red[i].outputs);
    }
    free(test);
    free(blue);
    free(red);
}

// find the inputs of the candidate components
void findInputs(int *sol_blue, int *sol_red)
{
    int i, j, k, aux, aux2, i_l, i_h, comp_size, n_reduc, *sol_blue_reduc,
        *sol_red_reduc, *sol_blue_reduc_t, *sol_red_reduc_t;
    gate_structure gate;

    for (i = 0; i < n_cand; i++) {
        comp_size = (int) ceil(size[i] / 2); // the inputs/outputs are created
                                             // with maximum size
        blue[i].inputs = new_gate_structure(comp_size);
        blue[i].outputs = new_gate_structure(comp_size);
        blue[i].first_entry.num = -1;
        red[i].inputs = new_gate_structure(comp_size);
        red[i].outputs = new_gate_structure(comp_size);
        red[i].first_entry.num = -1;
    }

    // Solutions without common edges
    sol_blue_reduc = alloc_vectori(gpx_n);
    sol_red_reduc = alloc_vectori(gpx_n);
    sol_blue_reduc_t = alloc_vectori(gpx_n);
    sol_red_reduc_t = alloc_vectori(gpx_n);
    j = k = 0;
    for (i = 0; i < gpx_n; i++) {
        aux = sol_blue[i];
        if (size[id[aux]] > 1) {
            sol_blue_reduc[j] = aux;
            sol_blue_reduc_t[j] = i;
            j++;
        }
        aux = sol_red[i];
        if (size[id[aux]] > 1) {
            sol_red_reduc[k] = aux;
            sol_red_reduc_t[k] = i;
            k++;
        }
    }
    n_reduc = j;

    // Blue Tour
    for (i = 0; i < n_cand; i++) {
        n_inputs[i] = 0;
        n_outputs[i] = 0;
        M_neigh[i][0] = 0; // the first column indicates the number
                           // of neighbours       
    }
    for (i = 0; i < n_reduc; i++) {
        gate.num = sol_blue_reduc[i];
        gate.time = sol_blue_reduc_t[i];
        if (i == 0)
            i_l = n_reduc - 1;
        else
            i_l = i - 1;
        if (i == (n_reduc - 1))
            i_h = 0;
        else
            i_h = i + 1;
        aux = id[sol_blue_reduc[i]];
        aux2 = id[sol_blue_reduc[i_l]];
        if (aux != aux2) {
            blue[aux].inputs[n_inputs[aux]] = gate;
            n_inputs[aux] = n_inputs[aux] + 1;
            if (blue[aux].first_entry.num == -1)
                blue[aux].first_entry = gate;
            // records the first entry in the candidate for the blue tour
            // Updating the neighbourhood relations
            if (M_neigh[aux][0] == 0) {
                M_neigh[aux][0] = 1;
                M_neigh[aux][1] = aux2;
                M_neigh2[aux][0] = 1;
            } else if (M_neigh[aux][0] == 1) {
                if (M_neigh[aux][1] == aux2)
                    M_neigh2[aux][0] = M_neigh2[aux][0] + 1;
                else {
                    M_neigh[aux][0] = 2;
                    M_neigh[aux][2] = aux2;
                    M_neigh2[aux][1] = 1;
                }
            } else {
                if (M_neigh[aux][1] == aux2)
                    M_neigh2[aux][0] = M_neigh2[aux][0] + 1;
                else if (M_neigh[aux][2] == aux2)
                    M_neigh2[aux][1] = M_neigh2[aux][1] + 1;
                else
                    M_neigh[aux][0] = M_neigh[aux][0] + 1;
            }
        }
        aux2 = id[sol_blue_reduc[i_h]];
        if (aux != aux2) {
            blue[aux].outputs[n_outputs[aux]] = gate;
            n_outputs[aux] = n_outputs[aux] + 1;
            blue[aux].last_exit = gate;
            // records the last exit in the candidate for the blue tour
            // Updating the neighbourhood relations
            if (M_neigh[aux][0] == 0) {
                M_neigh[aux][0] = 1;
                M_neigh[aux][1] = aux2;
                M_neigh2[aux][0] = 1;
            } else if (M_neigh[aux][0] == 1) {
                if (M_neigh[aux][1] == aux2)
                    M_neigh2[aux][0] = M_neigh2[aux][0] + 1;
                else {
                    M_neigh[aux][0] = 2;
                    M_neigh[aux][2] = aux2;
                    M_neigh2[aux][1] = 1;
                }
            } else {
                if (M_neigh[aux][1] == aux2)
                    M_neigh2[aux][0] = M_neigh2[aux][0] + 1;
                else if (M_neigh[aux][2] == aux2)
                    M_neigh2[aux][1] = M_neigh2[aux][1] + 1;
                else
                    M_neigh[aux][0] = M_neigh[aux][0] + 1;
            }
        }
    }

    // Red Tour
    for (i = 0; i < n_cand; i++) {
        n_inputs[i] = 0;
        n_outputs[i] = 0;
    }
    for (i = 0; i < n_reduc; i++) {
        gate.num = sol_red_reduc[i];
        gate.time = sol_red_reduc_t[i];
        if (i == 0)
            i_l = n_reduc - 1;
        else
            i_l = i - 1;
        if (i == (n_reduc - 1))
            i_h = 0;
        else
            i_h = i + 1;
        aux = id[sol_red_reduc[i]];
        if (aux != id[sol_red_reduc[i_l]]) {
            red[aux].inputs[n_inputs[aux]] = gate;
            n_inputs[aux] = n_inputs[aux] + 1;
            if (red[aux].first_entry.num == -1)
                red[aux].first_entry = gate;
            // records the first entry in the candidate for the red tour 
        }
        if (aux != id[sol_red_reduc[i_h]]) {
            red[aux].outputs[n_outputs[aux]] = gate;
            n_outputs[aux] = n_outputs[aux] + 1;
            red[aux].last_exit = gate;
            // records the last exit in the candidate for the red tour      
        }
    }

    free(sol_blue_reduc);
    free(sol_red_reduc);
    free(sol_blue_reduc_t);
    free(sol_red_reduc_t);
}

// test if graphs are equal (obs.: remember that each vertex has 0 or 1 edge)
int isequal(Graph * G1, Graph * G2)
{
    int i = 0, G1_empty, G2_empty, equal = 1;

    while (equal && i < G1->numVertices) {
        G1_empty = !G1->firstAdj[i]; // check if node i has an empty edge list
        G2_empty = !G2->firstAdj[i]; // check if node i has an empty edge list
        if (G1_empty != G2_empty)
            equal = 0;
        else {
            if (!G1_empty) {
                Adj *a = G1->firstAdj[i];
                // first edge of the adjacency list for node i
                Adj *b = G2->firstAdj[i];
                // first edge of the adjacency list for node i
                if (a->vertex != b->vertex)
                    equal = 0;
            }
        }
        i++;
    }
    return equal;
}

void testComp(int cand)
{
    int i, *inp_out_blue_inv, *inp_out_red;

    if (size[cand] <= 1)
        test[cand] = -1;
    else {
        if (n_inputs[cand] < 1)
            test[cand] = 0;
        else {
            if (n_inputs[cand] == 1)
                test[cand] = 1;
            else {
                // Graphs for blue and red tours                        
                inp_out_blue_inv = new_int(gpx_n);
                inp_out_red = new_int(2 * n_inputs[cand]);
                // Two cases: Case 1 - first input before first output;
                //            Case 2 - otherwise
                // Obs.: remember that the inputs and outputs are
                // inserted according to the flow                               
                if (blue[cand].inputs[0].time < blue[cand].outputs[0].time) {
                    for (i = 0; i < n_inputs[cand]; i++) {
                        inp_out_blue_inv[blue[cand].inputs[i].num] = 2 * i;
                        inp_out_blue_inv[blue[cand].outputs[i].num] =
                            2 * i + 1;
                    }
                } else {
                    for (i = 0; i < n_inputs[cand]; i++) {
                        inp_out_blue_inv[blue[cand].outputs[i].num] =
                            2 * i;
                        inp_out_blue_inv[blue[cand].inputs[i].num] =
                            2 * i + 1;
                    }
                }
                if (red[cand].inputs[0].time < red[cand].outputs[0].time) {
                    for (i = 0; i < n_inputs[cand]; i++) {
                        inp_out_red[2 * i] = red[cand].inputs[i].num;
                        inp_out_red[2 * i + 1] = red[cand].outputs[i].num;
                    }
                } else {
                    for (i = 0; i < n_inputs[cand]; i++) {
                        inp_out_red[2 * i] = red[cand].outputs[i].num;
                        inp_out_red[2 * i + 1] = red[cand].inputs[i].num;
                    }
                }

                Graph *Gs_blue = new_Graph(2 * n_inputs[cand]);
                // simplified graph for the blue path inside the component cand
                Graph *Gs_red = new_Graph(2 * n_inputs[cand]);
                // simplified graph for the red path inside the component cand
                // Two cases: Case 1 - first input before first output;
                //            Case 2 - otherwise
                // Obs.: remember that the inputs and outputs are inserted
                // according to the flow
                if (blue[cand].inputs[0].time < blue[cand].outputs[0].time) {
                    for (i = 0; i < (2 * n_inputs[cand] - 1); i = i + 2) {
                        insertEdge(Gs_blue, i, i + 1);
                        insertEdge(Gs_blue, i + 1, i);
                    }
                } else {
                    for (i = 1; i < (2 * n_inputs[cand] - 2); i = i + 2) {
                        insertEdge(Gs_blue, i, i + 1);
                        insertEdge(Gs_blue, i + 1, i);
                    }
                    insertEdge(Gs_blue, 2 * n_inputs[cand] - 1, 0);
                    insertEdge(Gs_blue, 0, 2 * n_inputs[cand] - 1);
                }
                if (red[cand].inputs[0].time < red[cand].outputs[0].time) {
                    for (i = 0; i < (2 * n_inputs[cand] - 1); i = i + 2) {
                        insertEdge(Gs_red,
                                   inp_out_blue_inv[inp_out_red[i]],
                                   inp_out_blue_inv[inp_out_red[i + 1]]);
                        insertEdge(Gs_red,
                                   inp_out_blue_inv[inp_out_red[i + 1]],
                                   inp_out_blue_inv[inp_out_red[i]]);
                    }
                } else {
                    for (i = 1; i < (2 * n_inputs[cand] - 2); i = i + 2) {
                        insertEdge(Gs_red,
                                   inp_out_blue_inv[inp_out_red[i]],
                                   inp_out_blue_inv[inp_out_red[i + 1]]);
                        insertEdge(Gs_red,
                                   inp_out_blue_inv[inp_out_red[i + 1]],
                                   inp_out_blue_inv[inp_out_red[i]]);
                    }
                    insertEdge(Gs_red,
                               inp_out_blue_inv[inp_out_red
                                                [2 * n_inputs[cand] - 1]],
                               inp_out_blue_inv[inp_out_red[0]]);
                    insertEdge(Gs_red,
                               inp_out_blue_inv[inp_out_red[0]],
                               inp_out_blue_inv[inp_out_red
                                                [2 * n_inputs[cand] - 1]]);
                }

                // Comparing the two graphs                                     
                test[cand] = isequal(Gs_blue, Gs_red);
                free(Gs_red);
                free(Gs_blue);
                free(inp_out_blue_inv);
                free(inp_out_red);
            }
        }
    }
}

// test if simplified graphs outside unfesible candidate component are equal
// Observation: this is equivalent of testing if all entries for a component
// are grouped after removing the feasible components (identified according to
// testComp) of the list of candidate entries
int testUnfeasibleComp(int *sol_blue)
{
    int i, j, aux, *comp_seq, *inp_comp_seq, n_newpart = 0;

    comp_seq = alloc_vectori(gpx_n); // sequence of all entries in unfeasible
                                 // components in the order given by sol_blue
    inp_comp_seq = alloc_vectori(n_cand);
    // records the number of entries in
    // each component in comp_seq
    // creating comp_seq
    j = 0; // j is the effective size of comp_seq
    aux = id[sol_blue[0]];
    if (test[aux] == 0 && aux != id[sol_blue[gpx_n - 1]])
        comp_seq[j++] = aux;
    inp_comp_seq[aux] = 0;
    for (i = 1; i < gpx_n; i++) {
        aux = id[sol_blue[i]];
        if (test[aux] == 0 && aux != id[sol_blue[i - 1]])
            comp_seq[j++] = aux;
        inp_comp_seq[aux] = 0;
    }

    // testing by checking the grouping of the components
    // (i.e., testing if the number of entries is 2
    if (j > 0) {
        aux = comp_seq[0];
        if (aux != comp_seq[j - 1])
            inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        for (i = 1; i < j; i++) {
            aux = comp_seq[i];
            if (aux != comp_seq[i - 1])
                inp_comp_seq[aux] = inp_comp_seq[aux] + 1;
        }
        for (i = 0; i < n_cand; i++)
            if (test[i] == 0 && inp_comp_seq[i] == 1) {
                test[i] = 1;
                n_newpart++;
            }
    }
    free(inp_comp_seq);
    free(comp_seq);
    return n_newpart;
}

// if candidate cand did not pass the test and has conditions,
// apply fusion with the neighbour with more connections
// (for more than 2 cutting points)
void fusion(int *sol_blue, int *sol_red)
{
    int i, cand, aux, n_fusions = 0, *neigh_vec_cond, *neigh_vec_ind;

    neigh_vec_cond = new_int(n_cand);
    // cond= -1 if cand is a true component or if it is between two
    // common edges;   = 0 stil does not chosen;
    //                 = 1 already chosen, and id will be changed;
    //                 = 2 already chosen but id will not be changed 
    neigh_vec_ind = new_int(n_cand); // neigh_vec_ind: indicates the neighbour
                                     // that will be fusioned with cand
    for (cand = 0; cand < n_cand; cand++)
        neigh_vec_cond[cand] = test[cand] == 1 || size[cand] <= 1 ? -1 : 0;

    for (cand = 0; cand < n_cand; cand++) {
        if (neigh_vec_cond[cand] == 0) {
            if (M_neigh[cand][0] == 1) {
                aux = M_neigh[cand][1];
                if (neigh_vec_cond[aux] == 0) {
                    neigh_vec_ind[cand] = aux;
                    neigh_vec_cond[cand] = 1;
                    neigh_vec_cond[aux] = 2;
                    n_fusions++;
                }
            } else if (M_neigh[cand][0] == 2) {
                if (M_neigh2[cand][0] > M_neigh2[cand][1]) {
                    aux = M_neigh[cand][1];
                    if (neigh_vec_cond[aux] == 0) {
                        neigh_vec_ind[cand] = aux;
                        neigh_vec_cond[cand] = 1;
                        neigh_vec_cond[aux] = 2;
                        n_fusions++;
                    }
                } else {
                    aux = M_neigh[cand][2];
                    if (neigh_vec_cond[aux] == 0) {
                        neigh_vec_ind[cand] = aux;
                        neigh_vec_cond[cand] = 1;
                        neigh_vec_cond[aux] = 2;
                        n_fusions++;
                    }
                }
            }

        }
    }

    if (n_fusions > 0) {
        // Reseting tour structures blue and red 
        for (cand = 0; cand < n_cand; cand++) {
            free(blue[cand].inputs);
            free(blue[cand].outputs);
            free(red[cand].inputs);
            free(red[cand].outputs);
        }
        free(blue);
        free(red);
        blue = new_tour(n_cand);
        red = new_tour(n_cand);

        // Making the fusions
        for (i = 0; i < gpx_n; i++) {
            aux = id[i];
            if (neigh_vec_cond[aux] == 1) {
                size[aux] = size[aux] - 1;
                id[i] = neigh_vec_ind[aux];
                aux = id[i];
                size[aux] = size[aux] + 1;
            }
        }

        // Repeating Step 5: Finding the inputs and outputs of each
        // candidate component
        findInputs(sol_blue, sol_red);
        // Repeating Step 6: testing the candidate components
        for (cand = 0; cand < n_cand; cand++)
            if (neigh_vec_cond[cand] == 2)
                testComp(cand); // test component cand                  
    }

    free(neigh_vec_cond);
    free(neigh_vec_ind);
}

// fusions of the candidate components in order to create partitions with two
// cutting points
void fusionB(int *sol_blue, int *sol_red)
{
    int i, cand, n_cand_seq = 0, previous_cand, next_cand, n_cand_new,
        n_rounds = 0, n_rounds_max = 1000;
    int *cand_seq, *cand_seq_cut, *assigned_cand, *new_label,
        *vector_new_cand, *new_component, n_newpart;

    // Memory allocation
    cand_seq = alloc_vectori(gpx_n); // list of entries and exits of
                                 // unfeasible candidates
    cand_seq_cut = alloc_vectori(gpx_n);
    new_label = alloc_vectori(n_cand);
    new_component = alloc_vectori(n_cand);
    Graph *G_cand = new_Graph(n_cand);

    // Walking in the blue tour and finding the high level cuts
    previous_cand = id[sol_blue[gpx_n - 1]];
    for (i = 0; i < gpx_n; i++) {
        cand = id[sol_blue[i]]; // candidate for vertex i of the blue tour
        if (i == (gpx_n - 1))
            next_cand = id[sol_blue[0]];
        else
            next_cand = id[sol_blue[i + 1]];
        // test if it is an unfeasible partition
        if (test[cand] == 0 && n_inputs[cand] > 0 && size[cand] > 1) {
            if (cand != previous_cand || cand != next_cand) {
                cand_seq[n_cand_seq] = cand;
                // checking if it is a first common entry or last common exit
                if (sol_blue[i] == blue[cand].first_entry.num &&
                    (blue[cand].first_entry.num ==
                     red[cand].first_entry.num
                     || blue[cand].first_entry.num ==
                     red[cand].last_exit.num))
                    cand_seq_cut[n_cand_seq] = 1; // it is a first common entry
                else if (sol_blue[i] == blue[cand].last_exit.num &&
                         (blue[cand].last_exit.num ==
                          red[cand].last_exit.num ||
                          blue[cand].last_exit.num ==
                          red[cand].first_entry.num))
                    cand_seq_cut[n_cand_seq] = 1; // it is a last common exit       
                else
                    cand_seq_cut[n_cand_seq] = 0;
                n_cand_seq++;
            }
        }
        previous_cand = cand;
    }
    // building the graph with conections between unfeasible components,
    // but without main entries and main exits
    if (n_cand_seq > 0) {
        for (i = 0; i < n_cand_seq - 1; i++) {
            if (cand_seq[i] != cand_seq[i + 1] &&
                cand_seq_cut[i] == 0 && cand_seq_cut[i + 1] == 0) {
                // insert edge between candidates
                insertEdge(G_cand, cand_seq[i], cand_seq[i + 1]);
                // insert edge between candidates                               
                insertEdge(G_cand, cand_seq[i + 1], cand_seq[i]);
            }
        }
        if (cand_seq[n_cand_seq - 1] != cand_seq[0] &&
            cand_seq_cut[n_cand_seq - 1] == 0 && cand_seq_cut[0] == 0) {
            // insert edge between candidates
            insertEdge(G_cand, cand_seq[n_cand_seq - 1], cand_seq[0]);
            // insert edge between candidates                           
            insertEdge(G_cand, cand_seq[0], cand_seq[n_cand_seq - 1]);
        }
    }
    // Walking in the red tour and finding the high level cuts
    n_cand_seq = 0;
    previous_cand = id[sol_red[gpx_n - 1]];
    for (i = 0; i < gpx_n; i++) {
        cand = id[sol_red[i]]; // candidate for vertex i of the red tour
        if (i == (gpx_n - 1))
            next_cand = id[sol_red[0]];
        else
            next_cand = id[sol_red[i + 1]];
        // test if it is an unfeasible partition
        if (test[cand] == 0 && n_inputs[cand] > 0 && size[cand] > 1) {
            if (cand != previous_cand || cand != next_cand) {
                cand_seq[n_cand_seq] = cand;
                // checking if it is a first common entry or last common exit
                if (sol_red[i] == red[cand].first_entry.num &&
                    (red[cand].first_entry.num ==
                     blue[cand].first_entry.num
                     || red[cand].first_entry.num ==
                     blue[cand].last_exit.num))
                    cand_seq_cut[n_cand_seq] = 1; // it is a first common entry
                else if (sol_red[i] == red[cand].last_exit.num &&
                         (red[cand].last_exit.num ==
                          blue[cand].last_exit.num ||
                          red[cand].last_exit.num ==
                          blue[cand].first_entry.num))
                    cand_seq_cut[n_cand_seq] = 1; // it is a last common exit       
                else
                    cand_seq_cut[n_cand_seq] = 0;
                n_cand_seq++;
            }
        }
        previous_cand = cand;
    }
    // building the graph with connections between unfeasible components,i
    // but without main entries and main exits
    if (n_cand_seq > 0) {
        for (i = 0; i < n_cand_seq - 1; i++) {
            if (cand_seq[i] != cand_seq[i + 1] &&
                cand_seq_cut[i] == 0 && cand_seq_cut[i + 1] == 0) {
                // insert edge between candidates
                insertEdge(G_cand, cand_seq[i], cand_seq[i + 1]);
                // insert edge between candidates                               
                insertEdge(G_cand, cand_seq[i + 1], cand_seq[i]);
            }
        }
        if (cand_seq[n_cand_seq - 1] != cand_seq[0] &&
            cand_seq_cut[n_cand_seq - 1] == 0 && cand_seq_cut[0] == 0) {
            // insert edge between candidates
            insertEdge(G_cand, cand_seq[n_cand_seq - 1], cand_seq[0]);
            // insert edge between candidates
            insertEdge(G_cand, cand_seq[0], cand_seq[n_cand_seq - 1]);
        }
    }

    for (cand = 0; cand < n_cand; cand++)
        new_label[cand] = cand;

    if (n_cand_seq > 0) {
        vector_new_cand = alloc_vectori(n_cand); // fusion of candidates        
        compCon(G_cand, vector_new_cand); // find the connected components of
                                          // the graph
        // new label
        n_cand_new = -1;
        for (i = 0; i < n_cand; i++) {
            new_component[i] = 0;
            if (n_cand_new < vector_new_cand[i])
                n_cand_new = vector_new_cand[i];
        }
        if (n_cand_new > -1) {
            n_cand_new++;
            assigned_cand = alloc_vectori(n_cand_new);
            for (cand = 0; cand < n_cand_new; cand++)
                assigned_cand[cand] = -1;
            for (cand = 0; cand < n_cand; cand++) {
                if (test[cand] == 0 && n_inputs[cand] > 0
                    && size[cand] > 1) {
                    if (assigned_cand[vector_new_cand[cand]] == -1) {
                        assigned_cand[vector_new_cand[cand]] = cand;
                        new_component[cand] = 1;
                    }
                    new_label[cand] = assigned_cand[vector_new_cand[cand]];
                }
            }
            free(assigned_cand);

            // Reseting tour structures blue and red 
            for (cand = 0; cand < n_cand; cand++) {
                free(blue[cand].inputs);
                free(blue[cand].outputs);
                free(red[cand].inputs);
                free(red[cand].outputs);
            }
            free(blue);
            free(red);
            blue = new_tour(n_cand);
            red = new_tour(n_cand);

            // Making the fusions
            for (i = 0; i < gpx_n; i++) {
                cand = id[i];
                if (new_label[cand] != cand) {
                    size[cand] = size[cand] - 1;
                    id[i] = new_label[cand];
                    cand = id[i];
                    size[cand] = size[cand] + 1;
                }
            }

            // Repeating Step 5: Finding the inputs and outputs of each
            // candidate component
            findInputs(sol_blue, sol_red);

            // Repeating Step 6: testing the candidate components
            // this procedure is O(n)
            for (cand = 0; cand < n_cand; cand++)
                if (new_component[cand] == 1)
                    testComp(cand);

            // Testing unfeasible partitions 
            do {
                n_rounds++;
                n_newpart = testUnfeasibleComp(sol_blue);
            } while (n_newpart > 0 && n_rounds <= n_rounds_max);
        }

        free(vector_new_cand);
    }
    free(G_cand);
    free(cand_seq);
    free(cand_seq_cut);
    free(new_label);
    free(new_component);
}

// select between the blue and red paths for each component
GainType off_gen(int *sol_blue, int *sol_red, int *offspring,
                 int *label_list)
{
    int i, i_l, i_h, k, aux, aux2, *select_cand, select_rest,
        *offspring_p2, *sol_blue_index, *sol_red_index, v_aux;
    GainType blue_fitness_rest, red_fitness_rest, fitness;

    // Computing the Fitness of each subtour in each true component
    for (i = 0; i < n_cand; i++) {
        blue[i].fitness = 0;
        red[i].fitness = 0;
    }
    blue_fitness_rest = 0; // fitness of the parts that are not true components
                           // in the blue path
    red_fitness_rest = 0;  // fitness of the parts that are not true components
                           // in the red path
    for (i = 0; i < gpx_n; i++) {
        if (i < gpx_n - 1)
            i_h = i + 1;
        else
            i_h = 0;
        aux = sol_blue[i];
        aux2 = sol_blue[i_h];
        if (test[id[aux]] > 0 && id[aux] == id[aux2]) {
            if (label_list[aux] != label_list[aux2])
                blue[id[aux]].fitness = blue[id[aux]].fitness +
                    weight(label_list[aux], label_list[aux2]);
        } else {
            if (label_list[aux] != label_list[aux2])
                blue_fitness_rest = blue_fitness_rest +
                    weight(label_list[aux], label_list[aux2]);
        }
        aux = sol_red[i];
        aux2 = sol_red[i_h];
        if (test[id[aux]] > 0 && id[aux] == id[aux2]) {
            if (label_list[aux] != label_list[aux2])
                red[id[aux]].fitness = red[id[aux]].fitness +
                    weight(label_list[aux], label_list[aux2]);
        } else if (label_list[aux] != label_list[aux2])
            red_fitness_rest = red_fitness_rest +
                weight(label_list[aux], label_list[aux2]);
    }

    // Selecting the components
    if (blue_fitness_rest > red_fitness_rest)
        select_rest = 1; // 1 for red and 0 for blue
    else
        select_rest = 0;
    select_cand = new_int(gpx_n);

    for (i = 0; i < n_cand; i++) {
        if (test[i] < 1)
            select_cand[i] = select_rest;
        else {
            if (blue[i].fitness > red[i].fitness)
                select_cand[i] = 1;
            else
                select_cand[i] = 0;
        }
    }

    // Generating the offspring     
    // index for the solutions
    sol_blue_index = new_int(gpx_n);
    sol_red_index = new_int(gpx_n);
    for (i = 0; i < gpx_n; i++) {
        sol_blue_index[sol_blue[i]] = i;
        sol_red_index[sol_red[i]] = i;
    }
    // Offspring: Graph
    Graph *G_offspring = new_Graph(gpx_n); // graph for offspring
    for (i = 0; i < gpx_n; i++) {
        if (select_cand[id[i]] == 0) {
            if (sol_blue_index[i] == 0) {
                i_l = sol_blue[gpx_n - 1];
                i_h = sol_blue[sol_blue_index[i] + 1];
            } else if (sol_blue_index[i] == (gpx_n - 1)) {
                i_l = sol_blue[sol_blue_index[i] - 1];
                i_h = sol_blue[0];
            } else {
                i_l = sol_blue[sol_blue_index[i] - 1];
                i_h = sol_blue[sol_blue_index[i] + 1];
            }
        } else {
            if (sol_red_index[i] == 0) {
                i_l = sol_red[gpx_n - 1];
                i_h = sol_red[sol_red_index[i] + 1];
            } else if (sol_red_index[i] == (gpx_n - 1)) {
                i_l = sol_red[sol_red_index[i] - 1];
                i_h = sol_red[0];
            } else {
                i_l = sol_red[sol_red_index[i] - 1];
                i_h = sol_red[sol_red_index[i] + 1];
            }
        }
        insertEdge(G_offspring, i, i_l);
        insertEdge(G_offspring, i, i_h);
    }

    // Offspring: vector
    offspring_p2 = new_int(gpx_n);
    i = 0;
    offspring_p2[i] = 0;
    v_aux = 0;
    Adj *a = G_offspring->firstAdj[v_aux]; // first edge of the vertex 
    Adj *b = a->nextAdj; // next edge of the vertex 
    v_aux = a->vertex;
    i++;
    offspring_p2[i] = v_aux;

    while (v_aux > 0) {
        a = G_offspring->firstAdj[v_aux]; // first edge of the vertex
        b = a->nextAdj; // next edge of the vertex
        assert(a->vertex != b->vertex);
        if (a->vertex == offspring_p2[i - 1])
            v_aux = b->vertex;
        else
            v_aux = a->vertex;
        i++;
        if (v_aux > 0)
            offspring_p2[i] = v_aux;
    }

    freeGraph(G_offspring);

    // Fitness of the offspring
    fitness = 0;
    for (i = 0; i < n_cand; i++) {
        if (test[i] == 1) {
            if (select_cand[i] == 0)
                fitness = fitness + blue[i].fitness;
            else
                fitness = fitness + red[i].fitness;
        }
    }
    if (select_rest == 0)
        fitness = fitness + blue_fitness_rest;
    else
        fitness = fitness + red_fitness_rest;

    // removing the ghost nodes
    k = 0;
    for (i = 0; i < gpx_n; i++)
        if (offspring_p2[i] < n_cities)
            offspring[k++] = offspring_p2[i];

    free(sol_red_index);
    free(sol_blue_index);
    free(select_cand);
    free(offspring_p2);
    return fitness;
}

/******************************************************************************\
*			 		 Dynamic allocation	and deallocation                       *
\******************************************************************************/

int *alloc_vectori(int lines)
{
    int *vector;

    vector = (int *) malloc(lines * sizeof(int));
    if (!vector) {
        printf("Allocation Error\n");
        exit(EXIT_FAILURE);
    }
    return vector;
}

int **alloc_matrixi(int lines, int collumns)
{
    int i, **Matrix;

    Matrix = (int **) malloc(lines * sizeof(int *));
    for (i = 0; i < lines; i++)
        Matrix[i] = (int *) malloc(collumns * sizeof(int));
    if (!Matrix) {
        printf("Allocation Error\n");
        exit(EXIT_FAILURE);
    }
    return Matrix;
}

void dealloc_matrixi(int **Matrix, int lines)
{
    int i;
    for (i = 0; i < lines; i++)
        free(Matrix[i]);
    free(Matrix);
}

/******************************************************************************\
*                                  Graph                                       *
\******************************************************************************/

Graph *new_Graph(int n)
{
    Graph *g = (Graph *) malloc(sizeof(Graph));
    g->numVertices = n;
    g->firstAdj = (Adj **) calloc(n, sizeof(Adj *));
    g->lastAdj = (Adj **) calloc(n, sizeof(Adj *));
    return g;
}

void insertEdge(Graph * g, int v1, int v2)
{
    Adj *a = (Adj *) malloc(sizeof(Adj));
    a->vertex = v2;
    a->nextAdj = 0;
    if (!g->firstAdj[v1])
        g->firstAdj[v1] = g->lastAdj[v1] = a;
    else
        g->lastAdj[v1] = g->lastAdj[v1]->nextAdj = a;
}

void freeGraph(Graph * g)
{
    int v;
    for (v = 0; v < g->numVertices; v++) {
        Adj *a = g->firstAdj[v];
        while (a) {
            Adj *tmp = a;
            a = a->nextAdj;
            free(tmp);
        }
    }
    free(g->firstAdj);
    free(g->lastAdj);
    free(g);
}

const int white = 0, grey = 1, black = 2;

void visitaDfsCC(Graph * g, int u, int *color,
                 int *vector_comp, int components)
{
    vector_comp[u] = components;
    color[u] = grey;
    Adj *a = g->firstAdj[u];
    while (a) {
        int v = a->vertex;
        if (color[v] == white)
            visitaDfsCC(g, v, color, vector_comp, components);
        a = a->nextAdj;
    }
    color[u] = black;
}

void compCon(Graph * g, int *vector_comp)
{
    int components = 0, u, gpx_n = g->numVertices;
    int *color = (int *) malloc(gpx_n * sizeof(int));

    for (u = 0; u < gpx_n; u++)
        color[u] = white;
    for (u = 0; u < gpx_n; u++) {
        if (color[u] == white) {
            visitaDfsCC(g, u, color, vector_comp, components);
            components++;
        }
    }
}
//################################### gpx.c end ###################################

//################################### Hashing.c begin ###################################
#include "Hashing.h"

/*
 * The functions HashInitialize, HashInsert and HashSearch is used
 * to maintain a hash table of tours. 
 *
 * A hash function maps tours to locations in a hash table. Each time 
 * a tour improvement has been found, the hash table is consulted to 
 * see whether the new tour happens to be local optimum found earlier. 
 * If this is the case, fruitless checkout time is avoided. 
 */

/*
 * HashInitialize(T) empties the hash table T.  
 * Empty entries have Cost equal to MINUS_INFINITY. 
 */

void HashInitialize(HashTable * T)
{
    int i;

    for (i = 0; i < HashTableSize; i++) {
        T->Entry[i].Hash = UINT_MAX;
        T->Entry[i].Cost = MINUS_INFINITY;
    }
    T->Count = 0;
}

/*
 * HashInsert(T,H,Cost) inserts H and Cost (the cost of the tour) in 
 * the table T in a location given by the hash value H. 
 *
 * Collisions are handled by double hashing.
 *
 * The table size is fixed. If the load factor becomes greater than 
 * a specified maximum, MaxLoadFactor, no more insertions will be
 * made. However, if the table entry given by H has a cost greater 
 * than or equal Cost, then Cost of this entry replaces its pervious 
 * value.      
 */

void HashInsert(HashTable * T, unsigned Hash, GainType Cost)
{
    int i = Hash % HashTableSize;
    if (T->Count >= MaxLoadFactor * HashTableSize) {
        if (Cost > T->Entry[i].Cost)
            return;
    } else {
        int p = Hash % 97 + 1;
        while (T->Entry[i].Cost != MINUS_INFINITY)
            if ((i -= p) < 0)
                i += HashTableSize;
        T->Count++;
    }
    T->Entry[i].Hash = Hash;
    T->Entry[i].Cost = Cost;
}

/*
 * HashSearch(T,H,Cost) returns 1 if table T has an entry containing 
 * Cost and H. Otherwise, the function returns 0.
 */

int HashSearch(HashTable * T, unsigned Hash, GainType Cost)
{
    int i, p;

    i = Hash % HashTableSize;
    p = Hash % 97 + 1;
    while ((T->Entry[i].Hash != Hash || T->Entry[i].Cost != Cost)
           && T->Entry[i].Cost != MINUS_INFINITY)
        if ((i -= p) < 0)
            i += HashTableSize;
    return T->Entry[i].Hash == Hash;
}
//################################### Hashing.c end ###################################

//################################### Heap.c begin ###################################
#include "LKH.h"
#include "Heap.h"

/*
 * A binary heap is used to implement a priority queue. 
 *
 * A heap is useful in order to speed up the computations of minimum 
 * spanning trees. The elements of the heap are the nodes, and the
 * priorities (ranks) are their associated costs (their minimum distance 
 * to the current tree). 
 */

static int HeapCount;           /* Its current number of elements */
static int HeapCapacity;        /* Its capacity */

/*      
 * The MakeHeap function creates an empty heap. 
 */

void HeapMake(int Size)
{
    Heap = (Node **) malloc((Size + 1) * sizeof(Node *));
    HeapCapacity = Size;
    HeapCount = 0;
}

/*
 * The HeapSiftUp function is called when the rank of a node is decreased, 
 * or when a node is inserted into the heap.
 * The function moves the node forward in the heap (the foremost node
 * of the heap has the lowest rank).
 * When calling HeapSiftUp(N), node N must belong to the heap.              
 */

void HeapSiftUp(Node * N)
{
    int Loc = N->Loc, Parent = Loc / 2;

    while (Parent && N->Rank < Heap[Parent]->Rank) {
        Heap[Loc] = Heap[Parent];
        Heap[Loc]->Loc = Loc;
        Loc = Parent;
        Parent /= 2;
    }
    Heap[Loc] = N;
    N->Loc = Loc;
}

/*
 * The HeapSiftDown function is called by the Heapify and HeapDeleteMin 
 * functions. The function moves the node backwards in the heap 
 * (the foremost node of the heap has the lowest rank).
 * When calling HeapSiftDown(N), node N must belong to the heap.              
 */

void HeapSiftDown(Node * N)
{
    int Loc = N->Loc, Child;

    while (Loc <= HeapCount / 2) {
        Child = 2 * Loc;
        if (Child < HeapCount && Heap[Child + 1]->Rank < Heap[Child]->Rank)
            Child++;
        if (N->Rank <= Heap[Child]->Rank)
            break;
        Heap[Loc] = Heap[Child];
        Heap[Loc]->Loc = Loc;
        Loc = Child;
    }
    Heap[Loc] = N;
    N->Loc = Loc;
}

/*       
 * The HeapDeleteMin function deletes the foremost node from the heap. 
 * The function returns a pointer to the deleted node (0, if the heap
 * is empty).
 */

Node *HeapDeleteMin()
{
    Node *Remove;

    if (!HeapCount)
        return 0;
    Remove = Heap[1];
    Heap[1] = Heap[HeapCount--];
    Heap[1]->Loc = 1;
    HeapSiftDown(Heap[1]);
    Remove->Loc = 0;
    return Remove;
}

/*       
 * The HeapInsert function inserts a node N into the heap.
 * When calling HeapInsert(N), node N must not belong to the heap.
 */

void HeapInsert(Node * N)
{
    HeapLazyInsert(N);
    HeapSiftUp(N);
}

/*
 * The HeapDelete function deletes a node N from the heap.
 */

void HeapDelete(Node * N)
{
    int Loc = N->Loc;
    if (!Loc)
        return;
    Heap[Loc] = Heap[HeapCount--];
    Heap[Loc]->Loc = Loc;
    if (Heap[Loc]->Rank > N->Rank)
        HeapSiftDown(Heap[Loc]);
    else
        HeapSiftUp(Heap[Loc]);
    N->Loc = 0;
}

/*       
 * The HeapLazyInsert function inserts a node as the last node of the heap.
 * This may destroy the heap condition, but it can later be restored by 
 * calling the Heapify function.
 * When calling HeapLazyInsert(N), node N must not belong to the heap.
 */

void HeapLazyInsert(Node * N)
{
    assert(HeapCount < HeapCapacity);
    Heap[++HeapCount] = N;
    N->Loc = HeapCount;
}

/*       
 * The Heapify function constructs a heap from its nodes.
 */

void Heapify()
{
    int Loc;
    for (Loc = HeapCount / 2; Loc >= 1; Loc--)
        HeapSiftDown(Heap[Loc]);
}
//################################### Heap.c end ###################################

//################################### IsBackboneCandidate.c begin ###################################
#include "LKH.h"

/* 
 * The IsBackboneCandidate function is used to test if an edge, (ta,tb), 
 * belongs to the set of backbone candidate edges.
 *
 * If the edge is a candidate edge the function returns 1; otherwise 0.
 */

int IsBackboneCandidate(const Node * ta, const Node * tb)
{
    Candidate *Nta;

    for (Nta = ta->BackboneCandidateSet; Nta && Nta->To; Nta++)
        if (Nta->To == tb)
            return 1;
    return 0;
}
//################################### IsBackboneCandidate.c end ###################################

//################################### IsCandidate.c begin ###################################
#include "LKH.h"

/* 
 * The IsCandidate function is used to test if an edge, (ta,tb), 
 * belongs to the set of candidate edges.
 *
 * If the edge is a candidate edge the function returns 1; otherwise 0.
 */

int IsCandidate(const Node * ta, const Node * tb)
{
    Candidate *Nta;

    for (Nta = ta->CandidateSet; Nta && Nta->To; Nta++)
        if (Nta->To == tb)
            return 1;
    return 0;
}
//################################### IsCandidate.c end ###################################

//################################### IsCommonEdge.c begin ###################################
#include "LKH.h"

/* 
 * The IsCommonEdge function is used to test if an edge, (ta,tb), 
 * is common to the set of tours to be merged.
 *
 * If the edge is a common edge the function returns 1; otherwise 0.
 */

int IsCommonEdge(const Node * ta, const Node * tb)
{
    int i;

    if (MergeTourFiles < 2)
        return 0;
    for (i = 0; i < MergeTourFiles; i++)
        if (ta->MergeSuc[i] != tb && tb->MergeSuc[i] != ta)
            return 0;
    return 1;
}
//################################### IsCommonEdge.c end ###################################

//################################### IsPossibleCandidate.c begin ###################################
#include "LKH.h"

/*
 * The IsPossibleCandidate function is used to test if an edge, (From,To),
 * may be a solution edge together with all fixed or common edges.
 *
 * If the edge is possible, the function returns 1; otherwise 0.
 */

int IsPossibleCandidate(Node * From, Node * To)
{
    Node *Na, *Nb, *Nc, *N;

    if (Forbidden(From, To))
        return 0;
    if (InInitialTour(From, To) ||
        From->SubproblemSuc == To || To->SubproblemSuc == From ||
        FixedOrCommon(From, To))
        return 1;
    if (From->FixedTo2 || To->FixedTo2)
        return 0;
    if (!IsCandidate(From, To) &&
        (FixedOrCommonCandidates(From) == 2 ||
         FixedOrCommonCandidates(To) == 2))
        return 0;
    if (MergeTourFiles < 2)
        return 1;
    if (!From->Head) {
        Na = FirstNode;
        do
            Na->Head = Na->Tail = Na;
        while ((Na = Na->Suc) != FirstNode);
        while ((Nb = Na->MergeSuc[0]) != FirstNode
               && FixedOrCommon(Na, Nb))
            Na = Nb;
        if (Nb != FirstNode) {
            N = Nb;
            do {
                Nc = Nb;
                do {
                    Na = Nb;
                    Na->Head = Nc;
                    Nb = Na->MergeSuc[0];
                } while (FixedOrCommon(Na, Nb));
                do
                    Nc->Tail = Na;
                while (Nc != N && (Nc = Nc->MergeSuc[0]) != Nb);
            } while (Nc != N);
        } else {
            do
                Nb->Head = Nb->Tail = FirstNode;
            while ((Nb = Nb->Suc) != FirstNode);
        }
    }
    if (From->Head == To->Head ||
        (From->Head != From && From->Tail != From) ||
        (To->Head != To && To->Tail != To))
        return 0;
    return 1;
}
//################################### IsPossibleCandidate.c end ###################################

//################################### KSwapKick.c begin ###################################
#include "LKH.h"

/*
 * The KSwapKick function makes a random walk K-swap kick, K>=4
 * (an extension of the double-bridge kick).
 *
 * The algorithm is inspired by the thesis
 *
 *    D. Richter,
 *    Toleranzen in Helsgauns Lin-Kernighan.Heuristik fur das TSP,
 *    Diplomarbeit, Martin-Luther-Universitat Halle-Wittenberg, 2006.
 *
 * and the paper
 *
 *    D. Applegate, W. Cook, and A. Rohe,
 *    Chained Lin-Kernighan for Large Traveling Salesman Problems.
 *    INFORMS Journal on Computing, 15 (1), pp. 82-92, 2003.
 */

#define WALK_STEPS 100
#define HUNT_COUNT (10 + Dimension / 1000)

static Node *RandomNode();
static Node *RandomWalkNode(Node * N);
static Node *LongEdgeNode();
static int KSwapKick_compare(const void *Na, const void *Nb);

void KSwapKick(int K)
{
    Node **s, *N;
    int Count, i;

    s = (Node **) malloc(K * sizeof(Node *));
    Count = 0;
    N = FirstNode;
    do {
        N->Rank = ++Count;
        N->V = 0;
    } while ((N = N->Suc) != FirstNode);
    N = LongEdgeNode();
    if (!N)
        goto End_KSwapKick;
    FirstNode = s[0] = N;
    N->V = 1;
    for (i = 1; i < K; i++) {
        N = s[i] = RandomWalkNode(s[0]);
        if (!N)
            K = i;
        else
            N->V = 1;
    }
    if (K < 4)
        goto End_KSwapKick;
    qsort(s, K, sizeof(Node *), KSwapKick_compare);
    for (i = 0; i < K; i++)
        s[i]->OldSuc = s[i]->Suc;
    for (i = 0; i < K; i++)
        Link(s[(i + 2) % K], s[i]->OldSuc);
  End_KSwapKick:
    free(s);
}

/*
 * The RandomNode function returns a random node N, for
 * which the edge (N, N->Suc) is neither a fixed edge nor
 * a common edge of tours to be merged, and N has not
 * previously been chosen.
 */

static Node *RandomNode()
{
    Node *N;
    int Count;

    if (SubproblemSize == 0)
        N = &NodeSet[1 + Random() % Dimension];
    else {
        N = FirstNode;
        for (Count = Random() % Dimension; Count > 0; Count--)
            N = N->Suc;
    }
    Count = 0;
    while ((FixedOrCommon(N, N->Suc) || N->V) && Count < Dimension) {
        N = N->Suc;
        Count++;
    }
    return Count < Dimension ? N : 0;
}

/*
 * The RandomWalkNode function makes a random walk on the
 * candidate edges starting at node N and returns a random
 * node R, for which the edge (R, R->Suc) is neither a
 * fixed edge nor a common edge of tours to be merged,
 * and R has not previously been chosen.
 */

static Node *RandomWalkNode(Node * N)
{
    Node *R = 0, *Last = 0;
    Candidate *NN;
    int Count, i;

    for (i = 1; i <= WALK_STEPS; i++) {
        Count = 0;
        for (NN = N->CandidateSet; NN->To; NN++)
            Count++;
        Count = Random() % Count;
        for (NN = N->CandidateSet; --Count > 0; NN++);
        if (NN->To != Last) {
            Last = N;
            N = NN->To;
            if (!N->V && !FixedOrCommon(N, N->Suc))
               R = N;
        }
    }
    return R ? R : RandomNode();
}

/*
 * The LongEdgeNode function condiders a small fraction of the nodes
 * and returns the one, R, that maximizes
 *
 *        C(R, R->Suc) - C(R, Nearest(R))
 */

static Node *LongEdgeNode()
{
    Node *N, *R = 0;
    int MaxG = INT_MIN, G, i;

    for (i = HUNT_COUNT; i > 0 && (N = RandomNode()); i--) {
        if ((G = C(N, N->Suc) - N->Cost) > MaxG) {
            MaxG = G;
            R = N;
        }
    }
    return R ? R : RandomNode();
}

static int KSwapKick_compare(const void *Na, const void *Nb)
{
    return (*(Node **) Na)->Rank - (*(Node **) Nb)->Rank;
}

//################################### KSwapKick.c end ###################################

//################################### LinKernighan.c begin ###################################
#include "Segment.h"
#include "LKH.h"
#include "Hashing.h"
#include "Sequence.h"

/*
 * The LinKernighan function seeks to improve a tour by sequential
 * and non-sequential edge exchanges.
 *
 * The function returns the cost of the resulting tour.
 */

GainType LinKernighan()
{
    Node *t1, *t2, *SUCt1;
    GainType Gain, G0, Cost;
    int X2, i, it = 0;
    Candidate *Nt1;
    Segment *S;
    SSegment *SS;
    double EntryTime = GetTime();

    Reversed = 0;
    S = FirstSegment;
    i = 0;
    do {
        S->Size = 0;
        S->Rank = ++i;
        S->Reversed = 0;
        S->First = S->Last = 0;
    }
    while ((S = S->Suc) != FirstSegment);
    SS = FirstSSegment;
    i = 0;
    do {
        SS->Size = 0;
        SS->Rank = ++i;
        SS->Reversed = 0;
        SS->First = SS->Last = 0;
    }
    while ((SS = SS->Suc) != FirstSSegment);

    FirstActive = LastActive = 0;
    Swaps = 0;

    /* Compute the cost of the initial tour, Cost.
       Compute the corresponding hash value, Hash.
       Initialize the segment list.
       Make all nodes "active" (so that they can be used as t1). */
    Cost = 0;
    Hash = 0;
    i = 0;
    t1 = FirstNode;
    do {
        t2 = t1->OldSuc = t1->Suc;
        t1->OldPred = t1->Pred;
        t1->Rank = ++i;
        Cost += (t1->SucCost = t2->PredCost = C(t1, t2)) - t1->Pi - t2->Pi;
        Hash ^= Rand[t1->Id] * Rand[t2->Id];
        t1->Cost = INT_MAX;
        for (Nt1 = t1->CandidateSet; (t2 = Nt1->To); Nt1++)
            if (t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
                t1->Cost = Nt1->Cost;
        t1->Parent = S;
        S->Size++;
        if (S->Size == 1)
            S->First = t1;
        S->Last = t1;
        if (SS->Size == 0)
            SS->First = S;
        S->Parent = SS;
        SS->Last = S;
        if (S->Size == GroupSize) {
            S = S->Suc;
            SS->Size++;
            if (SS->Size == SGroupSize)
                SS = SS->Suc;
        }
        t1->OldPredExcluded = t1->OldSucExcluded = 0;
        t1->Next = 0;
        if (KickType == 0 || Kicks == 0 || Trial == 1 ||
            !InBestTour(t1, t1->Pred) || !InBestTour(t1, t1->Suc))
            Activate(t1);
    }
    while ((t1 = t1->Suc) != FirstNode);
    if (S->Size < GroupSize)
        SS->Size++;
    Cost /= Precision;
    if (TraceLevel >= 3 || (TraceLevel == 2 && Cost < BetterCost)) {
        printff("Cost = " GainFormat, Cost);
        if (Optimum != MINUS_INFINITY && Optimum != 0)
            printff(", Gap = %0.4f%%", 100.0 * (Cost - Optimum) / Optimum);
        printff(", Time = %0.2f sec. %s\n", fabs(GetTime() - EntryTime),
                Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
    }
    PredSucCostAvailable = 1;

    /* Loop as long as improvements are found */
    do {
        /* Choose t1 as the first "active" node */
        while ((t1 = RemoveFirstActive())) {
            // if (GetTime() - EntryTime >= TimeLimit ||
            //     GetTime() - StartTime >= TotalTimeLimit) {
            //     if (TraceLevel >= 1)
            //         printff("*** Time limit exceeded in LinKernighan\n");
            //     goto End_LinKernighan;
            // }
            /* t1 is now "passive" */
            SUCt1 = SUC(t1);
            if ((TraceLevel >= 3 || (TraceLevel == 2 && Trial == 1)) &&
                ++it % (Dimension >= 100000 ? 10000 :
                        Dimension >= 10000 ? 1000 : 100) == 0)
                printff("#%d: Time = %0.2f sec.\n",
                        it, fabs(GetTime() - EntryTime));
            /* Choose t2 as one of t1's two neighbors on the tour */
            for (X2 = 1; X2 <= 2; X2++) {
                t2 = X2 == 1 ? PRED(t1) : SUCt1;
                if (FixedOrCommon(t1, t2) ||
                    (RestrictedSearch && Near(t1, t2) &&
                     (Trial == 1 ||
                      (Trial > BackboneTrials &&
                       (KickType == 0 || Kicks == 0)))))
                    continue;
                G0 = C(t1, t2);
                Gain = 0;
                /* Try to find a tour-improving chain of moves */
                do
                    t2 = Swaps == 0 ? BestMove(t1, t2, &G0, &Gain) :
                        BestSubsequentMove(t1, t2, &G0, &Gain);
                while (t2);
                if (Gain > 0) {
                    /* An improvement has been found */
#ifdef HAVE_LONG_LONG
                    assert(Gain % Precision == 0);
#else
                    assert(fmod(Gain, Precision) == 0);
#endif
                    Cost -= Gain / Precision;
                    if (TraceLevel >= 3 ||
                        (TraceLevel == 2 && Cost < BetterCost)) {
                        printff("Cost = " GainFormat, Cost);
                        if (Optimum != MINUS_INFINITY && Optimum != 0)
                            printff(", Gap = %0.4f%%",
                                    100.0 * (Cost - Optimum) / Optimum);
                        printff(", Time = %0.2f sec. %s\n",
                                fabs(GetTime() - EntryTime),
                                Cost < Optimum ? "<" : Cost ==
                                Optimum ? "=" : "");
                    }
                    StoreTour();
                    if (HashSearch(HTable, Hash, Cost))
                        goto End_LinKernighan;
                    /* Make t1 "active" again */
                    Activate(t1);
                    break;
                }
                RestoreTour();
                if (Dimension != DimensionSaved && SUC(t1) != SUCt1)
                    Reversed ^= 1;
            }
        }
        if (HashSearch(HTable, Hash, Cost))
            goto End_LinKernighan;
        HashInsert(HTable, Hash, Cost);
        /* Try to find improvements using non-sequential 4/5-opt moves */
        Gain = 0;
        if (Gain23Used && (Gain = Gain23()) > 0) {
            /* An improvement has been found */
#ifdef HAVE_LONG_LONG
            assert(Gain % Precision == 0);
#else
            assert(fmod(Gain, Precision) == 0);
#endif
            Cost -= Gain / Precision;
            StoreTour();
            if (TraceLevel >= 3 || (TraceLevel == 2 && Cost < BetterCost)) {
                printff("Cost = " GainFormat, Cost);
                if (Optimum != MINUS_INFINITY && Optimum != 0)
                    printff(", Gap = %0.4f%%",
                            100.0 * (Cost - Optimum) / Optimum);
                printff(", Time = %0.2f sec. + %s\n",
                        fabs(GetTime() - EntryTime),
                        Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
            }
            if (HashSearch(HTable, Hash, Cost))
                goto End_LinKernighan;

            // 移到这里判断，保证 StoreTour 执行至少一次
            // if (GetTime() - EntryTime >= TimeLimit ||
            //     GetTime() - StartTime >= TotalTimeLimit) {
            //     if (TraceLevel >= 1)
            //         printff("*** Time limit exceeded in LinKernighan\n");
            //     goto End_LinKernighan;
            // }
        }
    }
    while (Gain > 0);

  End_LinKernighan:
    PredSucCostAvailable = 0;
    NormalizeNodeList();
    NormalizeSegmentList();
    return Cost;
}
//################################### LinKernighan.c end ###################################

//################################### LKH.c begin ###################################
#include "LKH.h"
#include "Genetic.h"
#include "Sequence.h"
#include "gpx.h"

/* All global variables of the program. */

/* LKH.h variables */
// 新增
double TimeSpan; /* 统计改进值的时间跨度 */
double ScheduleScoreInSecond; /* 调度加分/秒 (未乘10) */
double SubProblemTotalTimeLimit; /* 分配给每个子问题的求解时间 */
double SubProblemStartTime; /* 当前子问题的开始求解时间 */

int AscentCandidates;   /* Number of candidate edges to be associated
                           with each node during the ascent */
int BackboneTrials;     /* Number of backbone trials in each run */
int Backtracking;       /* Specifies whether backtracking is used for 
                           the first move in a sequence of moves */
GainType BestCost;      /* Cost of the tour in BestTour */
int *BestTour;  /* Table containing best tour found */
GainType BetterCost;    /* Cost of the tour stored in BetterTour */
int *BetterTour;        /* Table containing the currently best tour 
                           in a run */
int CacheMask;  /* Mask for indexing the cache */
int *CacheVal;  /* Table of cached distances */
int *CacheSig;  /* Table of the signatures of cached 
                   distances */
int CandidateFiles;     /* Number of CANDIDATE_FILEs */
int *CostMatrix;        /* Cost matrix */
int Dimension;  /* Number of nodes in the problem */
int DimensionSaved;     /* Saved value of Dimension */
int EdgeFiles;          /* Number of EDGE_FILEs */
double Excess;  /* Maximum alpha-value allowed for any 
                   candidate edge is set to Excess times the 
                   absolute value of the lower bound of a 
                   solution tour */
int ExtraCandidates;    /* Number of extra neighbors to be added to 
                           the candidate set of each node */
Node *FirstActive, *LastActive; /* First and last node in the list 
                                   of "active" nodes */
Node *FirstNode;        /* First node in the list of nodes */
Segment *FirstSegment;  /* A pointer to the first segment in the cyclic 
                           list of segments */
SSegment *FirstSSegment;        /* A pointer to the first super segment in
                                   the cyclic list of segments */
int Gain23Used; /* Specifies whether Gain23 is used */
int GainCriterionUsed;  /* Specifies whether L&K's gain criterion is 
                           used */
double GridSize;        /* The grid size of toroidal instances */
int GroupSize;  /* Desired initial size of each segment */
int SGroupSize; /* Desired initial size of each super segment */
int Groups;     /* Current number of segments */
int SGroups;    /* Current number of super segments */
unsigned Hash;  /* Hash value corresponding to the current tour */
Node **Heap;    /* Heap used for computing minimum spanning 
                   trees */
HashTable *HTable;      /* Hash table used for storing tours */
int InitialPeriod;      /* Length of the first period in the ascent */
int InitialStepSize;    /* Initial step size used in the ascent */
double InitialTourFraction;     /* Fraction of the initial tour to be 
                                   constructed by INITIAL_TOUR_FILE edges */
char *LastLine; /* Last input line */
double LowerBound;      /* Lower bound found by the ascent */
int Kicks;      /* Specifies the number of K-swap-kicks */
int KickType;   /* Specifies K for a K-swap-kick */
int M;          /* The M-value is used when solving an ATSP-
                   instance by transforming it to a STSP-instance */
int MaxBreadth; /* The maximum number of candidate edges 
                   considered at each level of the search for
                   a move */
int MaxCandidates;      /* Maximum number of candidate edges to be 
                           associated with each node */
int MaxMatrixDimension; /* Maximum dimension for an explicit cost matrix */
int MaxSwaps;   /* Maximum number of swaps made during the 
                   search for a move */
int MaxTrials;  /* Maximum number of trials in each run */
int MergeTourFiles;     /* Number of MERGE_TOUR_FILEs */
int MoveType;   /* Specifies the sequantial move type to be used 
                   in local search. A value K >= 2 signifies 
                   that a k-opt moves are tried for k <= K */
Node *NodeSet;  /* Array of all nodes */
int Norm;       /* Measure of a 1-tree's discrepancy from a tour */
int NonsequentialMoveType;      /* Specifies the nonsequential move type to
                                   be used in local search. A value 
                                   L >= 4 signifies that nonsequential
                                   l-opt moves are tried for l <= L */
GainType Optimum;       /* Known optimal tour length. 
                           If StopAtOptimum is 1, a run will be 
                           terminated as soon as a tour length 
                           becomes equal this value */
int PatchingA;  /* Specifies the maximum number of alternating
                   cycles to be used for patching disjunct cycles */
int PatchingC;  /* Specifies the maximum number of disjoint cycles to be 
                   patched (by one or more alternating cycles) */
int Precision;  /* Internal precision in the representation of 
                   transformed distances */
int PredSucCostAvailable;  /* PredCost and SucCost are available */
int POPMUSIC_InitialTour;  /* Specifies whether the first POPMUSIC tour
                              is used as initial tour for LK */
int POPMUSIC_MaxNeighbors; /* Maximum number of nearest neighbors used 
                              as candidates in iterated 3-opt */
int POPMUSIC_SampleSize;   /* The sample size */
int POPMUSIC_Solutions;    /* Number of solutions to generate */
int POPMUSIC_Trials;       /* Maximum trials used for iterated 3-opt */
unsigned *Rand; /* Table of random values */
int Recombination; /* IPT or GPX2 */
int RestrictedSearch;      /* Specifies whether the choice of the first 
                              edge to be broken is restricted */
short Reversed; /* Boolean used to indicate whether a tour has 
                   been reversed */
int Run;        /* Current run number */
int Runs;       /* Total number of runs */
unsigned Seed;  /* Initial seed for random number generation */
double StartTime;       /* Time when execution starts */
int StopAtOptimum;      /* Specifies whether a run will be terminated if 
                           the tour length becomes equal to Optimum */
int Subgradient;        /* Specifies whether the Pi-values should be 
                           determined by subgradient optimization */
int SubproblemSize;     /* Number of nodes in a subproblem */
int SubsequentMoveType; /* Specifies the move type to be used for all 
                           moves following the first move in a sequence 
                           of moves. The value K >= 2 signifies that a 
                           K-opt move is to be used */
int SubsequentPatching; /* Species whether patching is used for 
                           subsequent moves */
SwapRecord *SwapStack;  /* Stack of SwapRecords */
int Swaps;      /* Number of swaps made during a tentative move */
double TimeLimit;       /* Time limit in seconds for each run */
double TotalTimeLimit;  /* Total time limit in seconds */
int TraceLevel; /* Specifies the level of detail of the output 
                   given during the solution process. 
                   The value 0 signifies a minimum amount of 
                   output. The higher the value is the more 
                   information is given */
int Trial;      /* Ordinal number of the current trial */

/* The following variables are read by the functions ReadParameters and 
   ReadProblem: */

char *ParameterFileName, *ProblemFileName, *PiFileName,
    *TourFileName, *OutputTourFileName, *InputTourFileName,
    **CandidateFileName, **EdgeFileName, *InitialTourFileName,
    *SubproblemTourFileName, **MergeTourFileName;
char *Name, *Type, *EdgeWeightType, *EdgeWeightFormat,
    *EdgeDataFormat, *NodeCoordType, *DisplayDataType;
int CandidateSetSymmetric, CandidateSetType,
    CoordType, DelaunayPartitioning, DelaunayPure,
    ExtraCandidateSetSymmetric, ExtraCandidateSetType,
    InitialTourAlgorithm,
    KarpPartitioning, KCenterPartitioning, KMeansPartitioning,
    MoorePartitioning,
    PatchingAExtended, PatchingARestricted,
    PatchingCExtended, PatchingCRestricted,
    ProblemType,
    RohePartitioning, SierpinskiPartitioning,
    SubproblemBorders, SubproblemsCompressed, WeightType, WeightFormat;

FILE *ParameterFile, *ProblemFile, *PiFile, *InputTourFile,
    *TourFile, *InitialTourFile, *SubproblemTourFile, **MergeTourFile;
CostFunction Distance, D, C, c;
MoveFunction BestMove, BacktrackMove, BestSubsequentMove;
MergeTourFunction MergeWithTour;

/* Genetic.h variables */
int MaxPopulationSize; /* The maximum size of the population */
int PopulationSize;    /* The current size of the population */
CrossoverFunction Crossover;
int **Population;      /* Array of individuals (solution tours) */
GainType *PenaltyFitness;  /* The fitnessl (tour penalty) of 
                              each individual */
GainType *Fitness;     /* The fitness (tour cost) of each individual */

/* Sequence.h variables */
Node **t;      /* The sequence of nodes to be used in a move */
Node **T;      /* The currently best t's */
Node **tSaved; /* For saving t when using the BacktrackKOptMove function */
int *p;        /* The permutation corresponding to the sequence in which
                  the t's occur on the tour */
int *q;        /* The inverse permutation of p */
int *incl;     /* Array: incl[i] == j, if (t[i], t[j]) is an
                  inclusion edge */
int *cycle;    /* Array: cycle[i] is cycle number of t[i] */
GainType *G;   /* For storing the G-values in the BestKOptMove
                  function */
int K;         /* The value K for the current K-opt move */

/* gpx.h variables */
int n_cities, n_cand;
int n_partitions_size2, n_partitions_before_fusion,
    n_partitions_after_fusion1, n_partitions_after_fusion2,
    n_partitions_after_fusion3;
int n_partitions_after_fusion4, n_partitions_after_fusion5,
    n_partitions_after_fusionB;
Node **Map2Node;
//################################### LKH.c end ###################################

//################################### LKHmain.c begin ###################################
#include "LKH.h"
#include "Genetic.h"

/*
 * This file contains the main function of the program.
 */

int LKHmain(double startTime)
{
    GainType Cost, OldOptimum;
    double Time, LastTime;

    StartTime = startTime;
    LastTime = GetTime();

    AllocateStructures();
    CreateCandidateSet();
    InitializeStatistics();

    if (Norm != 0)
        BestCost = PLUS_INFINITY;
    else {
        /* The ascent has solved the problem! */
        Optimum = BestCost = (GainType) LowerBound;
        UpdateStatistics(Optimum, GetTime() - LastTime);
        RecordBetterTour();
        RecordBestTour();
        Runs = 0;
    }

    /* Find a specified number (Runs) of local optima */
    for (Run = 1; Run <= Runs; Run++) {
        LastTime = GetTime();
        // if (Run!=1 && LastTime - StartTime >= TotalTimeLimit) {
        //     if (TraceLevel >= 1)
        //         printff("*** Time limit exceeded in LKHmain ***\n");
        //     Run--;
        //     break;
        // }
        Cost = FindTour();      /* using the Lin-Kernighan heuristic */
        
        if (Cost < BestCost) {
            BestCost = Cost;
            RecordBetterTour();
            RecordBestTour();
        }
        OldOptimum = Optimum;
        if (Cost < Optimum) {
            if (FirstNode->InputSuc) {
                Node *N = FirstNode;
                while ((N = N->InputSuc = N->Suc) != FirstNode);
            }
            Optimum = Cost;
            printff("*** New optimum = " GainFormat " ***\n", Optimum);
        }
        Time = fabs(GetTime() - LastTime);
        UpdateStatistics(Cost, Time);
        if (TraceLevel >= 1 && Cost != PLUS_INFINITY) {
            printff("Run %d: Cost = " GainFormat, Run, Cost);
            if (Optimum != MINUS_INFINITY && Optimum != 0)
                printff(", Gap = %0.4f%%",
                        100.0 * (Cost - Optimum) / Optimum);
            printff(", Time = %0.2f sec. %s\n\n", Time,
                    Cost < Optimum ? "<" : Cost == Optimum ? "=" : "");
        }
        SRandom(++Seed);
    }
    PrintStatistics();
    return EXIT_SUCCESS;
}
//################################### LKHmain.c end ###################################

//################################### Make2OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Make2OptMove function makes a 2-opt move by calling the macro Swap1 
 * (i.e., by calling either Flip of Flip_SL). Edges (t1,t2) and (t3,t4) 
 * are exchanged with edges (t2,t3) and (t4,t1). Node t4 is one of t3's 
 * two neighbors on the tour; which one is uniquely determined by the 
 * orientation of (t1,t2).
 */

void Make2OptMove(Node * t1, Node * t2, Node * t3, Node * t4)
{
    Swap1(t1, t2, t3);
}
//################################### Make2OptMove.c end ###################################

//################################### Make3OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Make3OptMove function makes a 3-opt move by calling the macro Swap2 
 * or Swap3.
 */

void
Make3OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
             Node * t5, Node * t6, int Case)
{
    switch (Case) {
    case 1:
    case 2:
        Swap2(t1, t2, t3, t6, t5, t4);
        return;
    case 5:
        Swap3(t1, t2, t4, t6, t5, t4, t6, t2, t3);
        return;
    case 6:
        Swap2(t3, t4, t5, t1, t2, t3);
        return;
    default:
        eprintf("Make3OptMove: Internal error");
    }
}
//################################### Make3OptMove.c end ###################################

//################################### Make4OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Make4OptMove function makes a 4-opt move by calling the macro Swap3.
 */

void
Make4OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
             Node * t5, Node * t6, Node * t7, Node * t8, int Case)
{
    if (SUC(t1) != t2)
        Reversed ^= 1;
    switch (Case) {
    case 1:
    case 2:
        Swap3(t1, t2, t3, t6, t5, t4, t7, t8, t1);
        return;
    case 3:
    case 4:
        Swap3(t1, t2, t3, t8, t7, t6, t5, t8, t1);
        return;
    case 5:
        if (!BETWEEN(t2, t7, t3))
            Swap3(t5, t6, t7, t2, t1, t4, t1, t4, t5);
        else if (BETWEEN(t2, t7, t6))
            Swap3(t5, t6, t7, t5, t8, t3, t3, t8, t1);
        else
            Swap3(t1, t2, t7, t7, t2, t3, t4, t7, t6);
        return;
    case 6:
        Swap3(t3, t4, t5, t6, t3, t2, t1, t6, t7);
        return;
    case 7:
        Swap3(t6, t5, t8, t2, t1, t4, t8, t5, t4);
        return;
    case 11:
        Swap3(t1, t2, t7, t3, t4, t5, t3, t6, t7);
        return;
    case 12:
        Swap3(t3, t4, t5, t7, t8, t1, t3, t6, t7);
        return;
    case 15:
        Swap3(t3, t4, t5, t3, t6, t7, t8, t3, t2);
        return;
    default:
        eprintf("Make4OptMove: Internal error");
    }
}
//################################### Make4OptMove.c end ###################################

//################################### Make5OptMove.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The Make5OptMove function makes a 5-opt move by calling the macro Swap4 
 * or Swap5.
 */

void
Make5OptMove(Node * t1, Node * t2, Node * t3, Node * t4,
             Node * t5, Node * t6, Node * t7, Node * t8,
             Node * t9, Node * t10, int Case)
{
    if (SUC(t1) != t2)
        Reversed ^= 1;
    switch (Case) {
    case 1:
        Swap4(t1, t2, t3, t8, t7, t6, t10, t9, t8, t10, t5, t4);
        return;
    case 2:
        if (BETWEEN(t2, t9, t4))
            Swap4(t1, t2, t3, t5, t6, t7, t10, t9, t8, t5, t10, t1);
        else
            Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t7, t10, t1);
        return;
    case 3:
        Swap4(t3, t4, t5, t7, t8, t9, t1, t2, t3, t7, t10, t1);
        return;
    case 4:
        Swap5(t5, t6, t8, t1, t2, t3, t10, t9, t8, t1, t4, t5, t6, t10,
              t1);
        return;
    case 5:
        Swap5(t5, t6, t10, t1, t2, t3, t6, t10, t1, t8, t7, t6, t8, t4,
              t5);
        return;
    case 6:
        Swap4(t1, t2, t3, t9, t10, t1, t7, t8, t9, t6, t5, t4);
        return;
    case 7:
        if (BETWEEN(t3, t9, t7))
            Swap4(t3, t4, t5, t8, t7, t6, t10, t9, t8, t1, t2, t3);
        else if (BETWEEN(t6, t9, t4))
            Swap4(t3, t4, t5, t8, t7, t6, t9, t10, t1, t9, t2, t3);
        else
            Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t9, t7, t10, t1);
        return;
    case 8:
        Swap4(t3, t4, t5, t9, t10, t1, t8, t7, t6, t8, t3, t2);
        return;
    case 9:
        Swap4(t10, t9, t8, t5, t6, t7, t1, t2, t3, t1, t4, t5);
        return;
    case 10:
        if (BETWEEN(t5, t9, t7))
            Swap4(t5, t6, t7, t9, t10, t1, t4, t3, t2, t4, t9, t8);
        else if (BETWEEN(t3, t9, t6))
            Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t9, t7, t10, t1);
        else
            Swap4(t1, t2, t3, t9, t10, t1, t5, t6, t7, t5, t8, t9);
        return;
    case 11:
        if (BETWEEN(t3, t9, t6))
            Swap4(t1, t2, t3, t6, t5, t4, t9, t10, t1, t7, t8, t9);
        else
            Swap4(t5, t6, t7, t10, t9, t8, t2, t1, t10, t4, t3, t2);
        return;
    case 12:
        Swap4(t1, t2, t3, t8, t7, t6, t10, t9, t8, t5, t10, t1);
        return;
    case 13:
        if (BETWEEN(t4, t9, t7))
            Swap5(t7, t8, t10, t5, t6, t7, t1, t2, t3, t5, t9, t1, t9, t1,
                  t10);
        else if (BETWEEN(t6, t9, t3))
            Swap5(t10, t9, t1, t7, t8, t9, t3, t4, t5, t3, t6, t7, t3, t1,
                  t10);
        else
            Swap5(t10, t9, t1, t4, t3, t2, t5, t6, t7, t5, t8, t10, t9, t1,
                  t10);
        return;
    case 14:
        Swap5(t10, t9, t1, t5, t6, t7, t5, t8, t9, t3, t4, t5, t3, t1,
              t10);
        return;
    case 15:
        if (BETWEEN(t6, t9, t3))
            Swap5(t10, t9, t1, t3, t4, t5, t6, t3, t2, t8, t7, t6, t9, t1,
                  t10);
        else
            Swap5(t1, t2, t6, t3, t4, t5, t8, t7, t6, t10, t9, t8, t2, t10,
                  t1);
        return;
    case 16:
        if (BETWEEN(t4, t9, t7))
            Swap4(t3, t4, t5, t8, t7, t6, t9, t10, t1, t8, t3, t2);
        else if (BETWEEN(t5, t9, t3))
            Swap4(t3, t4, t5, t9, t10, t1, t6, t3, t2, t7, t8, t9);
        else
            Swap4(t3, t4, t5, t1, t2, t3, t7, t8, t9, t7, t10, t1);
        return;
    case 17:
        if (BETWEEN(t7, t9, t3))
            Swap4(t3, t4, t5, t7, t8, t9, t2, t1, t10, t3, t6, t7);
        else
            Swap4(t7, t8, t9, t2, t1, t10, t3, t4, t5, t3, t6, t7);
        return;
    case 18:
        Swap4(t3, t4, t5, t7, t8, t9, t3, t6, t7, t1, t2, t3);
        return;
    case 19:
        Swap4(t7, t8, t9, t1, t2, t3, t6, t5, t4, t7, t10, t1);
        return;
    case 20:
        Swap4(t7, t8, t9, t3, t4, t5, t10, t7, t6, t3, t10, t1);
        return;
    case 21:
        Swap4(t5, t6, t7, t5, t8, t9, t1, t2, t3, t4, t1, t10);
        return;
    case 22:
        Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t1, t9, t10, t1);
        return;
    case 23:
        Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t1, t9, t10, t1);
        return;
    case 24:
        Swap4(t1, t2, t3, t8, t7, t6, t5, t8, t1, t9, t10, t1);
        return;
    case 25:
        Swap4(t1, t2, t3, t8, t7, t6, t5, t8, t1, t9, t10, t1);
        return;
    case 26:
        if (!BETWEEN(t2, t7, t3))
            Swap4(t5, t6, t7, t2, t1, t4, t1, t4, t5, t9, t10, t1);
        else if (BETWEEN(t2, t7, t6))
            Swap4(t5, t6, t7, t5, t8, t3, t3, t8, t1, t9, t10, t1);
        else
            Swap4(t1, t2, t7, t7, t2, t3, t4, t7, t6, t9, t10, t1);
        return;
    case 27:
        Swap4(t3, t4, t5, t6, t3, t2, t1, t6, t7, t9, t10, t1);
        return;
    case 28:
        Swap4(t6, t5, t8, t2, t1, t4, t8, t5, t4, t9, t10, t1);
        return;
    case 29:
        Swap4(t1, t2, t7, t3, t4, t5, t3, t6, t7, t9, t10, t1);
        return;
    case 30:
        if (BETWEEN(t3, t7, t5))
            Swap4(t3, t4, t5, t7, t8, t1, t7, t2, t3, t9, t10, t1);
        else
            Swap4(t3, t4, t5, t3, t6, t7, t1, t2, t3, t9, t10, t1);
        return;
    case 31:
        Swap4(t3, t4, t5, t3, t6, t7, t8, t3, t2, t9, t10, t1);
        return;
    case 32:
        Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t7, t10, t1);
        return;
    case 33:
        if (BETWEEN(t3, t9, t5))
            Swap4(t1, t2, t3, t5, t6, t7, t10, t9, t8, t5, t10, t1);
        else
            Swap4(t1, t2, t3, t7, t8, t9, t7, t10, t1, t5, t6, t7);
        return;
    case 34:
        Swap4(t7, t8, t9, t1, t2, t3, t1, t4, t5, t7, t10, t1);
        return;
    case 35:
        Swap4(t9, t10, t1, t5, t6, t7, t4, t3, t2, t9, t4, t5);
        return;
    case 36:
        Swap4(t9, t10, t1, t7, t8, t9, t3, t4, t5, t6, t3, t2);
        return;
    case 37:
        if (BETWEEN(t6, t9, t4))
            Swap4(t1, t2, t3, t6, t5, t4, t9, t10, t1, t8, t7, t6);
        else
            Swap4(t9, t10, t1, t3, t4, t5, t3, t6, t7, t3, t8, t9);
        return;
    case 38:
        if (BETWEEN(t3, t9, t7))
            Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t6, t1, t10);
        else if (BETWEEN(t6, t9, t4))
            Swap4(t1, t2, t3, t6, t5, t4, t7, t8, t9, t7, t10, t1);
        else
            Swap4(t3, t4, t5, t9, t10, t1, t8, t7, t6, t3, t8, t9);
        return;
    case 39:
        Swap4(t1, t2, t3, t7, t8, t9, t5, t6, t7, t1, t4, t5);
        return;
    case 40:
        Swap4(t9, t10, t1, t4, t3, t2, t5, t6, t7, t5, t8, t9);
        return;
    case 41:
        if (BETWEEN(t5, t9, t7))
            Swap4(t7, t8, t9, t1, t2, t3, t6, t5, t4, t7, t10, t1);
        else if (BETWEEN(t3, t9, t6))
            Swap4(t1, t2, t3, t5, t6, t7, t9, t10, t1, t5, t8, t9);
        else
            Swap4(t5, t6, t7, t9, t10, t1, t2, t9, t8, t3, t4, t5);
        return;
    case 42:
        if (BETWEEN(t3, t9, t6))
            Swap4(t7, t8, t9, t5, t6, t7, t1, t2, t3, t1, t4, t5);
        else
            Swap4(t9, t10, t1, t5, t6, t7, t3, t4, t5, t3, t8, t9);
        return;
    case 43:
        Swap4(t1, t2, t3, t7, t8, t9, t6, t5, t4, t7, t10, t1);
        return;
    case 44:
        if (BETWEEN(t4, t9, t7))
            Swap4(t7, t8, t9, t5, t6, t7, t1, t2, t3, t5, t10, t1);
        else if (BETWEEN(t6, t9, t3))
            Swap4(t9, t10, t1, t5, t6, t7, t3, t4, t5, t3, t8, t9);
        else
            Swap4(t7, t8, t9, t1, t2, t3, t6, t5, t4, t7, t10, t1);
        return;
    case 45:
        Swap4(t9, t10, t1, t3, t4, t5, t7, t8, t9, t3, t6, t7);
        return;
    case 46:
        Swap4(t7, t8, t9, t5, t6, t7, t3, t4, t5, t1, t2, t3);
        return;
    case 47:
        if (BETWEEN(t4, t9, t7))
            Swap4(t5, t6, t7, t1, t2, t3, t9, t10, t1, t5, t8, t9);
        else if (BETWEEN(t5, t9, t3))
            Swap4(t9, t10, t1, t7, t8, t9, t5, t6, t7, t3, t4, t5);
        else
            Swap4(t7, t8, t9, t3, t4, t5, t3, t6, t7, t2, t1, t10);
        return;
    case 48:
        if (BETWEEN(t7, t9, t3))
            Swap4(t3, t4, t5, t8, t7, t6, t2, t1, t10, t8, t3, t2);
        else
            Swap4(t3, t4, t5, t7, t8, t9, t3, t6, t7, t1, t2, t3);
        return;
    case 49:
        Swap4(t9, t10, t1, t5, t6, t7, t3, t4, t5, t3, t8, t9);
        return;
    case 50:
        Swap4(t3, t4, t5, t3, t6, t7, t9, t10, t1, t8, t3, t2);
        return;
    case 51:
        Swap4(t5, t6, t7, t1, t2, t3, t9, t10, t1, t4, t9, t8);
        return;
    case 52:
        Swap4(t5, t6, t7, t3, t4, t5, t9, t10, t1, t3, t8, t9);
        return;
    default:
        eprintf("Make5OptMove: Internal error");
    }
}
//################################### Make5OptMove.c end ###################################

//################################### MakeKOptMove.c begin ###################################
#include "Sequence.h"
#include "Segment.h"

/*
 * The MakeKOptMove function makes a K-opt move (K >= 2) using sorting by 
 * reversals.
 *   
 * Let t[1:2K] be the sequence of nodes used in the K-opt move.  
 * 
 *    {(t[2i-1],t[2i]) | 1 <= i <= K} is the set of edges to be excluded 
 *                                    from the tour. 
 
 *    {(t[2i],t[incl[2i]]) | 1 <= i <= K} is the set of edges to be 
 *                                        included in the tour.
 *   
 * And let p[1:2K] be a permutation corresponding to the sequence in which 
 * the nodes occur on the tour.
 *   
 * Then the move corresponds to sorting p by reversals. 
 *   
 * MakeKOptMove finds the minimum number of reversals and makes the 
 * corresponding series of 2-opt moves (swaps).
 *   
 * The implementation is based upon the algorithm for sorting signed 
 * permutations by reversals given in
 *   
 *    A, Bergeron,
 *    "A Very Elementary Presentation of the Hannenhalli-Pevzner Theory",
 *    Lecture Notes in Computer Science, 2089, 106-117 (2001). 
 */

static void Reverse(int i, int j);
static int Score(int Left, int Right, int K);

void MakeKOptMove(int K)
{
    int i, j, Best_i = 0, Best_j = 0, BestScore, s;

    FindPermutation(K);
  FindNextReversal:
    /* Find the oriented reversal that has maximal score */
    BestScore = -1;
    for (i = 1; i <= 2 * K - 2; i++) {
        j = q[incl[p[i]]];
        if (j >= i + 2 && (i & 1) == (j & 1) &&
            (s = i & 1 ? Score(i + 1, j, K) :
             Score(i, j - 1, K)) > BestScore) {
            BestScore = s;
            Best_i = i;
            Best_j = j;
        }
    }
    if (BestScore >= 0) {
        i = Best_i;
        j = Best_j;
        if (i & 1) {
            Swap1(t[p[i + 1]], t[p[i]], t[p[j]]);
            Reverse(i + 1, j);
        } else {
            Swap1(t[p[i - 1]], t[p[i]], t[p[j]]);
            Reverse(i, j - 1);
        }
        goto FindNextReversal;
    }
    /* No more oriented reversals. Cut a simpe hurdle, if any.
     * Note that there can be no super hurdles */
    for (i = 1; i <= 2 * K - 3; i += 2) {
        j = q[incl[p[i]]];
        if (j >= i + 3) {
            Swap1(t[p[i]], t[p[i + 1]], t[p[j]]);
            Reverse(i + 1, j - 1);
            goto FindNextReversal;
        }
    }
}

/*
 * The Reverse function reverses the sequence of elements in p[i:j].
 * The inverse permutation q is updated accordingly.
 */

static void Reverse(int i, int j)
{
    for (; i < j; i++, j--) {
        int pi = p[i];
        q[p[i] = p[j]] = i;
        q[p[j] = pi] = j;
    }
}

/*
 * The Score function computes the score of a reversal. The score is the 
 * number of oriented pairs in the resulting reversal.
 */

static int Score(int Left, int Right, int K)
{
    int Count = 0, i, j;

    Reverse(Left, Right);
    for (i = 1; i <= 2 * K - 2; i++) {
        j = q[incl[p[i]]];
        if (j >= i + 2 && (i & 1) == (j & 1))
            Count++;
    }
    Reverse(Left, Right);
    return Count;
}
//################################### MakeKOptMove.c end ###################################

//################################### MergeTourWithBestTour.c begin ###################################
#include "LKH.h"

/*
 * The MergeTourWithBestTour function attempts to find a short 
 * tour by merging the current tour with the tour in the array BestTour.
 * 
 * If a tour shorter than BestTour is found, Pred and Suc of each 
 * node point to its neighbors, and the tour cost is returned.
 */

GainType MergeTourWithBestTour()
{
    Node *N1, *N2, *M1, *M2;
    int i;
    // ATSP
    int Dim = Dimension / 2;
    for (i = 1; i <= Dim; i++) {
        N1 = &NodeSet[BestTour[i - 1]];
        N2 = &NodeSet[BestTour[i]];
        M1 = &NodeSet[N1->Id + Dim];
        M2 = &NodeSet[N2->Id + Dim];
        M1->Next = N1;
        N1->Next = M2;
        M2->Next = N2;
    }
    return MergeWithTour();
}
//################################### MergeTourWithBestTour.c end ###################################

//################################### MergeWithTourCLARIST.c begin ###################################
#include "LKH.h"
#include "CLARIST.h"

/*
 * The MergeWithTourCLARIST function attempts to find a short tour by merging
 * a given tour, T1, with another tour, T2. T1 is given by the Suc pointers
 * of its nodes. T2 is given by the Next pointers of its nodes.
 *
 * Originally programmed by X. Clarist. Adapted for LKH by K. Helsgaun.
 *
 * The merging algorithm may be described as follows: 
 * Let G be the graph consisting of the nodes and the union of the edges of
 * T1 and T2.
 * The algorithm first identifies the nodes of degree 2 in G.
 * It subsequently identifies the components. A component consists of
 * adjacent nodes in G that are not of degree 2. 
 * This base being established, the goal is to determine if a path of T1
 * that is part of a component can be replaced by a path of T2 that is part
 * of the same component, in the case this reduces the length of T1.
 * It may happen that such a replacement does not result in a tour,
 * because sometimes two or more components are related to each other,
 * and the replacement will only result in a tour if those components are
 * replaced all at once. Therefore, the algorithm tries to find the
 * neighboring components in order to fuse them. 
 * The following process is applied to determine if two components C1 and C2
 * may be fused:
 * Departing from node N of C1, let's move along the paths of T2 inside of
 * C1 and along the paths of T1 outside of C1. If N is met again, and only
 * C1 and C2 have been visited, they are fused. After a round of fusion on
 * components, the following validation process is used to determine if the
 * replacement of the paths of T1 by the paths of T2 inside of a component
 * (or a fused component) will result in a tour:
 * Let n1 be the number of accesses to the component, if the paths of T1
 * are followed inside the component and outside of it. 
 * Let n2 be the number of accesses to the component, if the paths of T2
 * are followed inside the component and the paths of T1 outside of it.
 * If n1 = n2, the replacement will result in a tour.
 * Iterations of fusions and validations are performed as many times as
 * possible (the valid components are excluded after each iteration).
 * Finally, the replacement of valid components that reduce the length is
 * made in T1. 
 *
 * Sometimes a shorter tour can be found by switching T1 and T2.
 * The algorithm above is therefore executed again after a switch of
 * T1 and T2, and the shortest tour found is selected.
 */

GainType MergeWithTourCLARIST()
{
    Node *N, *Prev;
    rec *ptcur;
    int len, i;
    GainType Cost1 = 0, Cost2 = 0, OldCost1, OldCost2;
    static GainType BestCost = PLUS_INFINITY;

    if (vecpttra == NULL) {
        int Dim = Dimension > DimensionSaved ? Dimension : DimensionSaved;
        vecpttra = (rec *) malloc((Dim + 1) * sizeof(rec));
        for (i = 1; i <= Dim; i++)
            vecpttra[i].ID = i;
        lnkdif = (double *) malloc((MAXDIFNBR + 1) * sizeof(double));
        lnkgrp = (double *) malloc((MAXDIFNBR + 1) * sizeof(double));
        grp2 = (int *) malloc((MAXDIFNBR + 1) * sizeof(int));
        grp2N = (int *) malloc((MAXDIFNBR + 1) * sizeof(int));
        difact = (int *) malloc((MAXDIFNBR + 2) * sizeof(int));
        difact++;
        diftst1 = (int *) malloc((MAXDIFNBR + 1) * sizeof(int));
        diftst2 = (int *) malloc((MAXDIFNBR + 1) * sizeof(int));
    }

    N = FirstNode;
    ptdeb = &vecpttra[N->Id];
    do {
        ptcur = &vecpttra[N->Id];
        (ptcur->ptN = ptcur->ptbufN = &vecpttra[N->Suc->Id])->ptP = ptcur;
        (ptcur->pt2N = &vecpttra[N->Next->Id])->pt2P = ptcur;
        ptcur->ptbufN->ptbufP = ptcur;
        N->Next->Prev = N->Suc->Pred = N;
        ptcur->len = (C(N, N->Suc) - N->Pi - N->Suc->Pi) / Precision;
        ptcur->len2 = (C(N, N->Next) - N->Pi - N->Next->Pi) / Precision;
        Cost1 += ptcur->len;
        Cost2 += ptcur->len2;
    } while ((N = N->Suc) != FirstNode);

    if (Cost1 == Cost2) {
        N = FirstNode;
        do {
            if (N->Suc != N->Next && N->Suc != N->Prev)
                break;
        } while ((N = N->Suc) != FirstNode);
        if (N == FirstNode &&
            (N->Suc == N->Next || N->Suc == N->Prev))
            return Cost1;
    }
    N = FirstNode;
    do
        N->OldSuc = N->Suc;
    while ((N = N->Suc) != FirstNode);
    OldCost1 = Cost1;
    OldCost2 = Cost2;

    Cost1 += merge_clarist();
    N = FirstNode;
    do {
        ptcur = &vecpttra[N->Id];
        ptcur->ptbufPsaved = 
            Cost1 <= Cost2 || Cost1 < OldCost1 ? ptcur->ptbufP
                                               : &vecpttra[N->Prev->Id];
        ptcur->ptbufNsaved =
            Cost1 <= Cost2 || Cost1 < OldCost1 ? ptcur->ptbufN
                                               : &vecpttra[N->Next->Id];
    } while ((N = N->Suc) != FirstNode);
    do {
        rec *pttmp;
        ptcur = &vecpttra[N->Id];
        pttmp = ptcur->pt2N;
        (ptcur->pt2N = ptcur->ptN)->pt2P = ptcur;
        (ptcur->ptN = ptcur->ptbufN = pttmp)->ptP = ptcur;
        ptcur->ptbufN->ptbufP = ptcur;
        len = ptcur->len;
        ptcur->len = ptcur->len2;
        ptcur->len2 = len;
    } while ((N = N->Suc) != FirstNode);
    BestCost = Cost1 <= Cost2 ? Cost1 : Cost2;
    Cost2 += merge_clarist();
    if (Cost2 < BestCost) {
        Prev = NULL;
        do {
            ptcur = &vecpttra[N->Id];
            N->Suc = &NodeSet[ptcur->ptbufN->ID] != Prev ?
                     &NodeSet[ptcur->ptbufN->ID] :
                     &NodeSet[ptcur->ptbufP->ID];
            N->Suc->Pred = Prev = N;
        } while ((N = N->Suc) != FirstNode);
        BestCost = Cost2;
    } else if (Cost1 == OldCost1 || Cost1 == OldCost2) {
        N = FirstNode;
        do
            (N->Suc = N->OldSuc)->Pred = N;
        while ((N = N->Suc) != FirstNode);
        return OldCost1;
    } else {
        Prev = NULL;
        do {
            ptcur = &vecpttra[N->Id];
            N->Suc = &NodeSet[ptcur->ptbufNsaved->ID] != Prev ?
                     &NodeSet[ptcur->ptbufNsaved->ID] :
                     &NodeSet[ptcur->ptbufPsaved->ID];
            N->Suc->Pred = Prev = N;
        } while ((N = N->Suc) != FirstNode);
        BestCost = Cost1;
    }
    Hash = 0;
    N = FirstNode;
    do
        Hash ^= Rand[N->Id] * Rand[N->Suc->Id];
    while ((N = N->Suc) != FirstNode);
    if (TraceLevel >= 2 && BestCost < OldCost1 && BestCost < OldCost2)
        printff("CLARIST: " GainFormat "\n", BestCost);
    return BestCost;
}

int merge_clarist()
{
    rec *ptcur;
    int i, j;

    reduce_path_tour1();
    tag_all_components();
    if (difnegfnd) {
        reduce_path_tour2();
        for (i = 1; i <= difnbr; i++) {
            grp2[i] = i;
            grp2N[i] = i;
            diftst1[i] = 0;
            lnkgrp[i] = lnkdif[i];
            difact[i] = 0;
        }
        j = 0;
        do {
            fusgrp2 = 0;
            j++;
            ptdebtog = ptdebcom2;
            fuse_components();
            if (j == 1 || fusgrp2)
                validate_components();
            stop = 1;
            if (fusgrp2) {
                for (i = 1; i <= difnbr; i++) {
                    if (!diftst1[grp2[i]]) {
                        stop = 0;
                        break;
                    }
                }
            }
            if (!stop) {
                ptcur = ptdebcom2;
                do
                    ptcur = ptcur->pt21->ptCC;
                while ((ptcur->pt21 == ptcur->pt22 ||
                        diftst1[grp2[ptcur->diftag]]) &&
                       ptcur != ptdebcom2);
                ptdebcom2 = ptcur;
                if (ptcur->pt21 != ptcur->pt22 &&
                    !diftst1[grp2[ptcur->diftag]]) {
                    do {
                        if (ptcur->pt21 == ptcur->pt22 ||
                            diftst1[grp2[ptcur->diftag]]) {
                            pttmp = ptcur->pt21->ptCC;
                            ptcur->ptCC->ptCC = pttmp;
                            pttmp->ptCC = ptcur->ptCC;
                            ptcur = pttmp;
                        } else
                            ptcur = ptcur->pt21->ptCC;
                    } while (ptcur != ptdebcom2);
                }
            }
        } while (!stop);
        totdif = 0;
        for (i = 1; i <= difnbr; i++) {
            if (diftst1[grp2[i]] && lnkgrp[grp2[i]] < 0) {
                difact[i] = 1;
                totdif += lnkdif[i];
            }
        }
        if (totdif < 0) {
            if (valid_tour()) {
                generate_offspring();
                return totdif;
            }
        }
    }
    return 0;
}

void reduce_path_tour1()
{
    rec *ptcur, *ptcom;

    ptcur = ptdeb;
    if (ptcur->pt2N != ptcur->ptN && ptcur->pt2P != ptcur->ptN) {
        do
            ptcur = ptcur->ptN;
        while (ptcur->pt2N != ptcur->ptN && ptcur->pt2P != ptcur->ptN);
    } else {
        while (ptcur->pt2N == ptcur->ptP || ptcur->pt2P == ptcur->ptP)
            ptcur = ptcur->ptP;
    }
    ptdebcom = ptcur;
    do {
        ptcom = ptcur;
        ptcur = ptcur->ptN;
        while (ptcur->pt2N == ptcur->ptN || ptcur->pt2P == ptcur->ptN) {
            ptcur->diftag = -1;
            ptcur = ptcur->ptN;
        }
        ptcur->ptC = ptcom;
        ptcom->ptC = ptcur;
        ptcur->diftag = 0;
        ptcur = ptcur->ptN;
        while (ptcur->pt2N != ptcur->ptN && ptcur->pt2P != ptcur->ptN) {
            ptcur->ptC = NULL;
            ptcur->diftag = 0;
            ptcur = ptcur->ptN;
        }
        ptcur->diftag = 0;
        ptcur->pt1C = ptcom;
    } while (ptcur != ptdebcom);
}

void find_component_extent(rec * ptcur)
{
    ptcur->diftag = difcnt;
    if (ptcur->ptP->diftag != -1)
        lnkcnt1 += ptcur->ptP->len;
    if (ptcur->ptN->diftag != -1)
        lnkcnt1 += ptcur->len;
    if (ptcur->pt2P->diftag != -1)
        lnkcnt2 += ptcur->pt2P->len2;
    if (ptcur->pt2N->diftag != -1)
        lnkcnt2 += ptcur->len2;
    if (ptcur->ptP->diftag == 0)
        find_component_extent(ptcur->ptP);
    if (ptcur->ptN->diftag == 0)
        find_component_extent(ptcur->ptN);
    if (ptcur->pt2P->diftag == 0)
        find_component_extent(ptcur->pt2P);
    if (ptcur->pt2N->diftag == 0)
        find_component_extent(ptcur->pt2N);
}

void tag_one_component(rec * ptcur)
{
    lnkcnt1 = 0;
    lnkcnt2 = 0;
    find_component_extent(ptcur);
    lnkdif[difcnt] = (lnkcnt2 - lnkcnt1) * 0.5;
    if (lnkdif[difcnt] < 0)
        difnegfnd = 1;
    difcnt++;
}

void tag_all_components()
{
    rec *ptcur, *ptcom, *pttmp;
    long long diftag;

    difcnt = 1;
    difnegfnd = 0;
    ptcur = ptdebcom;
    tag_one_component(ptcur);
    do
        ptcur = ptcur->pt1C;
    while (ptcur->diftag != 0 && ptcur != ptdebcom);
    totC21 = 0;
    ptdebcom2 = ptcur;

    if (ptcur->diftag == 0) {
        do {
            if (ptcur->diftag == 0)
                tag_one_component(ptcur);
            ptcom = ptcur;
            diftag = ptcom->diftag;
            do
                ptcur = ptcur->pt1C;
            while (ptcur->diftag == diftag);
            totC21++;
            pttmp = ptcur->ptC;
            pttmp->pt21 = ptcom;
            ptcom->pt21 = pttmp;
            ptcur->ptCC = pttmp;
            pttmp->ptCC = ptcur;
        } while (ptcur != ptdebcom2);
    } else {
        totC21++;
        pttmp = ptdebcom->ptC;
        pttmp->pt21 = ptdebcom;
        ptdebcom->pt21 = pttmp;
        ptdebcom->ptCC = pttmp;
        pttmp->ptCC = ptdebcom;
    }
    difnbr = difcnt - 1;
}

void reduce_path_tour2()
{
    rec *ptcur, *ptcom, *pttmp;
    long long diftag;

    ptcur = ptdebcom2;
    if (ptcur->ptN == ptcur->pt2N || ptcur->ptP == ptcur->pt2N)
        ptcur = ptcur->ptC;
    ptdebcom2 = ptcur;
    do {
        ptcom = ptcur;
        diftag = ptcom->diftag;
        do {
            do
                ptcur = ptcur->pt2N;
            while (ptcur->ptC == NULL);
            ptcur = ptcur->ptC;
        } while (ptcur->diftag == diftag && ptcur != ptdebcom2);
        pttmp = ptcur->ptC;
        pttmp->pt22 = ptcom;
        ptcom->pt22 = pttmp;
    } while (ptcur != ptdebcom2);
}

void fuse_components()
{
    rec *ptcur, *ptcom, *ptcom1;
    int diftag, idxdif, grp2idxdif;
    int unique;

    ptcom = ptcom1 = ptdebtog;
    do {
        idxdif = ptcom1->diftag;
        grp2idxdif = grp2[idxdif];
        ptcur = ptcom1->pt22->ptCC;
        diftag = 0;
        unique = 0;
        do {
            if (grp2[ptcur->diftag] != grp2idxdif) {
                if (grp2[ptcur->diftag] != grp2[diftag]) {
                    if (diftag == 0) {
                        diftag = ptcur->diftag;
                        unique = 1;
                    } else
                        unique = 0;
                }
                ptcur = ptcur->pt21->ptCC;
            } else
                ptcur = ptcur->pt22->ptCC;
        } while ((unique || diftag == 0) && ptcur != ptcom1);
        if (unique && diftag != 0) {
            fusgrp2 = 1;
            lnkgrp[grp2idxdif] += lnkgrp[grp2[diftag]];
            idx2 = diftag;
            grp2[idx2] = grp2idxdif;
            while (grp2N[idx2] != diftag) {
                idx2 = grp2N[idx2];
                grp2[idx2] = grp2idxdif;
            }
            grp2N[idx2] = grp2N[idxdif];
            grp2N[idxdif] = diftag;
        }
        do
            ptcom1 = ptcom1->pt21->ptCC;
        while (grp2[ptcom1->diftag] == grp2idxdif && ptcom1 != ptcom);
    } while (ptcom1 != ptcom);
}

void validate_components()
{
    rec *ptcom, *ptcur;
    rec *ptout[MAXDIFNBR + 1];
    int difcnt[MAXDIFNBR + 1];
    long long grp2diftag;
    int cnt, i;

    ptcur = ptdebcom2;
    do {
        ptcur->ptE = NULL;
        ptcur = ptcur->pt21;
        ptcur->ptE = NULL;
        ptcur = ptcur->ptCC;
    } while (ptcur != ptdebcom2);

    for (i = 1; i <= difnbr; i++) {
        ptout[i] = NULL;
        difcnt[i] = 0;
        diftst2[i] = 0;
    }
    ptcur = ptdebcom2;
    for (i = 1; i <= 2; i++) {
        do {
            grp2diftag = grp2[ptcur->diftag];
            pttmp = ptout[grp2diftag];
            if (pttmp != NULL) {
                ptcur->ptE = pttmp;
                pttmp->ptE = ptcur;
            }
            if (i == 2)
                difcnt[grp2diftag]++;
            ptcur = ptcur->pt21;
            ptout[grp2[ptcur->diftag]] = ptcur;
            ptcur = ptcur->ptCC;
        } while (ptcur != ptdebcom2);
    }
    ptcom = ptdebcom2;
    do {
        if (diftst2[grp2[ptcom->diftag]] == 0) {
            cnt = 0;
            ptcur = ptcom;
            do {
                cnt++;
                ptcur = ptcur->pt22->ptE;
            } while (ptcur != ptcom);
            grp2diftag = grp2[ptcom->diftag];
            diftst2[grp2diftag] = 1;
            if (cnt == difcnt[grp2diftag])
                diftst1[grp2diftag] = 1;
        }
        ptcom = ptcom->pt21->ptCC;
    } while (ptcom != ptdebcom2);
}

int valid_tour()
{
    rec *ptcur;

    cntC2 = 0;
    ptcur = ptdebcom2;
    do {
        if (!difact[ptcur->diftag])
            ptcur = ptcur->pt21->ptC;
        else
            ptcur = ptcur->pt22->ptC;
        cntC2++;
    } while (ptcur != ptdebcom2);
    return cntC2 == totC21;
}

void generate_offspring()
{
    rec *ptcur, *ptlas, *ptlas2, *pt1C;

    if (ptdebcom2->ptP == ptdebcom2->pt2P ||
        ptdebcom2->ptP == ptdebcom2->pt2N)
        ptdebcom2 = ptdebcom2->ptC;
    ptcur = ptdebcom2;
    do {
        if (difact[ptcur->diftag]) {
            ptlas2 = ptcur->pt21->ptC;
            do {
                pt1C = ptcur->pt1C;
                ptlas = pt1C->ptC->ptP;
                do {
                    ptcur->ptbufP = ptcur->pt2P;
                    ptcur->ptbufN = ptcur->pt2N;
                    ptcur->lenbuf = ptcur->len2;
                    ptcur = ptcur->ptP;
                } while (ptcur != ptlas);
                ptcur = pt1C;
            } while (ptcur != ptlas2);
        } else
            ptcur = ptcur->pt21->ptC;
    } while (ptcur != ptdebcom2);
}
//################################### MergeWithTourCLARIST.c end ###################################

//################################### MergeWithTourGPX2.c begin ###################################
#include "LKH.h"
#include "gpx.h"

/*
 * The MergeWithTourGPX2 function attempts to find a short tour
 * by merging a given tour, T1, with another tour, T2.
 * T1 is given by the Suc pointers of its nodes.
 * T2 is given by the Next pointers of its nodes.
 *
 * The merging algorithm uses Generalized Partition Crossover 2,
 * GPX2, described in
 *
 *      R.Tinos, D. Whitley, and G. Ochoa (2017),
 *      A new generalized partition crossover for the traveling
 *      salesman problem: tunneling between local optima.
 */

GainType MergeWithTourGPX2()
{
    int NewDimension = 0;
    GainType Cost1 = 0, ShrunkCost1 = 0, ShrunkCost2 = 0, NewCost;
    Node *N, *First = 0, *Last;
    int *red, *blue, *offspring, i;

    N = FirstNode;
    do
        N->Suc->Pred = N->Next->Prev = N;
    while ((N = N->Suc) != FirstNode);
    i = 0;
    do {
        Cost1 += C(N, N->Suc) - N->Pi - N->Suc->Pi;
        if ((N->Suc == N->Prev || N->Suc == N->Next) &&
            (N->Pred == N->Prev || N->Pred == N->Next))
            N->V = 0;
        else {
            N->V = 1;
            NewDimension++;
            First = N;
        }
    } while ((N = N->Suc) != FirstNode);
    Cost1 /= Precision;
    if (NewDimension == 0)
        return Cost1;

    /* Shrink the tours. 
       OldPred and OldSuc represent the shrunken T1. 
       Prev and Next represent the shrunken T2 */
    N = First;
    Last = 0;
    do {
        if (N->V) {
            if (Last)
                (Last->OldSuc = N)->OldPred = Last;
            Last = N;
        }
    } while ((N = N->Suc) != First);
    (Last->OldSuc = First)->OldPred = Last;
    Last = 0;
    do {
        if (N->V) {
            if (Last) {
                Last->Next = N;
                if (Last != N->Prev)
                    N->Prev = Last;
            }
            Last = N;
        }
    } while ((N = N->Next) != First);
    Last->Next = First;
    if (Last != First->Prev)
        First->Prev = Last;

    n_cities = NewDimension;
    red = (int *) malloc(n_cities * sizeof(int));
    blue = (int *) malloc(n_cities * sizeof(int));
    offspring = (int *) malloc((n_cities + 1) * sizeof(int));
    Map2Node = (Node **) malloc(n_cities * sizeof(Node *));

    N = First;
    i = 0;
    do {
        Map2Node[i] = N;
        red[i] = N->Rank = i;
        i++;
        ShrunkCost1 += C(N, N->OldSuc) - N->Pi - N->OldSuc->Pi;
    } while ((N = N->OldSuc) != First);
    i = 0;
    do {
        blue[i++] = N->Rank;
        ShrunkCost2 += C(N, N->Next) - N->Pi - N->Next->Pi;
    } while ((N = N->Next) != First);
    ShrunkCost1 /= Precision;
    ShrunkCost2 /= Precision;

    /* Perform GPX2 recombination */
    NewCost = gpx(red, blue, offspring);
    
    free(red);
    free(blue);
    if (NewCost >= ShrunkCost1 || NewCost >= ShrunkCost2) {
        free(offspring);
        free(Map2Node);
        return Cost1;
    }
    offspring[n_cities] = offspring[0];
    for (i = 0; i < n_cities; i++) {
        N = Map2Node[offspring[i]];
        Node *NextN = Map2Node[offspring[i + 1]];
        N->OldSuc = NextN;
        NextN->OldPred = N;
    }
    free(offspring);
    free(Map2Node);

    /* Expand the offspring into a full tour */
    N = FirstNode;
    do
        N->Mark = 0;
    while ((N = N->Suc) != FirstNode);
    N = First;
    N->Mark = N;
    do {
        if (!N->Suc->Mark && (!N->V || !N->Suc->V))
            N->OldSuc = N->Suc;
        else if (!N->Pred->Mark && (!N->V || !N->Pred->V))
            N->OldSuc = N->Pred;
        else if (N->OldSuc->Mark)
            N->OldSuc = !N->OldPred->Mark ? N->OldPred : First;
        N->Mark = N;
    } while ((N = N->OldSuc) != First);
    
    Cost1 = 0;
    Hash = 0;
    do {
        Cost1 += C(N, N->OldSuc) - N->Pi - N->OldSuc->Pi;
        N->OldSuc->Pred = N;
        Hash ^= Rand[N->Id] * Rand[N->OldSuc->Id];
    }
    while ((N = N->Suc = N->OldSuc) != First);
    Cost1 /= Precision;
    if (TraceLevel >= 2)
        printff("GPX2: " GainFormat "\n", Cost1);
    return Cost1;
}
//################################### MergeWithTourGPX2.c end ###################################

//################################### MergeWithTourIPT.c begin ###################################
#include "LKH.h"

/*
 * The MergeWithTourIPT function attempts to find a short tour
 * by merging a given tour, T1, with another tour, T2. 
 * T1 is given by the Suc pointers of its nodes. 
 * T2 is given by the Next pointers of its nodes.
 *
 * The merging algorithm may be described as follows:
 * Let G be the graph consisting of the nodes and the union of the 
 * edges of T1 and T2. Attempt - in all possible ways - 
 * to separate the nodes of G into two disjoint sets (A,B) such 
 * that the cardinality of edges connecting A with B is exactly two. 
 * If this is possible, any replacement of T1's A-edges with T2's
 * A-edges results in a tour. The same holds for the B-edges. 
 * If such a replacement reduces the cost of one the tours, then make it. 
 *
 * If a tour T shorter than T1 and T2 is found, Pred and Suc of each node
 * point to its neighbors in T, and T's cost is returned.
 *
 * The function is called from the FindTour function.
 *
 * The implementation is inspired by the algorithm described in the
 * paper
 *
 *   A. Mobius, B. Freisleben, P. Merz, and M. Schreiber, 
 *   "Combinatorial Optimization by Iterative Partial Transcription", 
 *   Physical Review E, Volume 59, Number 4, pp. 4667-4674, 1999.
 */

GainType MergeWithTourIPT()
{
    int Rank = 0, Improved1 = 0, Improved2 = 0;
    int SubSize1, SubSize2, MaxSubSize1, NewDimension = 0, Forward;
    int MinSubSize, BestMinSubSize = 3, MinForward = 0;
    GainType Cost1 = 0, Cost2 = 0, Gain, OldCost1, MinGain = 0;
    Node *N, *NNext, *N1, *N2, *MinN1, *MinN2, *First = 0, *Last;

    N = FirstNode;
    do
        N->Suc->Pred = N->Next->Prev = N;
    while ((N = N->Suc) != FirstNode);
    do {
        Cost1 += N->Cost = C(N, N->Suc) - N->Pi - N->Suc->Pi;
        if ((N->Suc == N->Prev || N->Suc == N->Next) &&
            (N->Pred == N->Prev || N->Pred == N->Next))
            N->V = 0;
        else {
            N->V = 1;
            NewDimension++;
            First = N;
        }
    } while ((N = N->Suc) != FirstNode);
    if (NewDimension == 0)
        return Cost1 / Precision;
    do {
        Cost2 += N->NextCost = N->Next == N->Pred ? N->Pred->Cost :
            N->Next == N->Suc ? N->Cost :
            C(N, N->Next) - N->Pi - N->Next->Pi;
    } while ((N = N->Next) != FirstNode);
    OldCost1 = Cost1;

    /* Shrink the tours. 
       OldPred and OldSuc represent the shrunken T1. 
       Prev and Next represent the shrunken T2 */
    N = First;
    Last = 0;
    do {
        if (N->V) {
            N->Rank = ++Rank;
            if (Last) {
                (Last->OldSuc = N)->OldPred = Last;
                if (Last != N->Pred)
                    Last->Cost = 0;
            }
            Last = N;
        }
    } while ((N = N->Suc) != First);
    (Last->OldSuc = First)->OldPred = Last;
    if (Last != First->Pred)
        Last->Cost = 0;
    N = First;
    Last = 0;
    do {
        if (N->V) {
            if (Last) {
                Last->Next = N;
                if (Last != N->Prev) {
                    N->Prev = Last;
                    Last->NextCost = 0;
                }
            }
            Last = N;
        }
    } while ((N = N->Next) != First);
    Last->Next = First;
    if (Last != First->Prev) {
        First->Prev = Last;
        Last->NextCost = 0;
    }

    /* Merge the shrunken tours */
    do {
        MinN1 = MinN2 = 0;
        MinSubSize = NewDimension / 2;
        N1 = First;
        do {
            while (N1->OldSuc != First &&
                   (N1->OldSuc == N1->Next || N1->OldSuc == N1->Prev))
                N1 = N1->OldSuc;
            if (N1->OldSuc == First &&
                (N1->OldSuc == N1->Next || N1->OldSuc == N1->Prev))
                break;
            for (Forward = 1, N2 = N1->Next; Forward >= 0;
                 Forward--, N2 = N1->Prev) {
                if (N2 == N1->OldSuc || N2 == N1->OldPred)
                    continue;
                SubSize2 = MaxSubSize1 = 0;
                do {
                    if (++SubSize2 >= MinSubSize)
                        break;
                    if ((SubSize1 = N2->Rank - N1->Rank) < 0)
                        SubSize1 += NewDimension;
                    if (SubSize1 >= MinSubSize)
                        break;
                    if (SubSize1 > MaxSubSize1) {
                        if (SubSize1 == SubSize2) {
                            for (N = N1, Gain = 0; N != N2; N = N->OldSuc)
                                Gain += N->Cost - N->NextCost;
                            if (!Forward)
                                Gain += N1->NextCost - N2->NextCost;
                            if (Gain != 0) {
                                MinSubSize = SubSize1;
                                MinN1 = N1;
                                MinN2 = N2;
                                MinGain = Gain;
                                MinForward = Forward;
                            }
                            break;
                        }
                        MaxSubSize1 = SubSize1;
                    }
                } while ((N2 = Forward ? N2->Next : N2->Prev) != N1);
            }
        } while ((N1 = N1->OldSuc) != First &&
                 MinSubSize != BestMinSubSize);
        if (MinN1) {
            BestMinSubSize = MinSubSize;
            if (MinGain > 0) {
                Improved1 = 1;
                Cost1 -= MinGain;
                Rank = MinN1->Rank;
                for (N = MinN1; N != MinN2; N = NNext) {
                    NNext = MinForward ? N->Next : N->Prev;
                    (N->OldSuc = NNext)->OldPred = N;
                    N->Rank = Rank;
                    N->Cost = MinForward ? N->NextCost : NNext->NextCost;
                    if (++Rank > NewDimension)
                        Rank = 1;
                }
            } else {
                Improved2 = 1;
                Cost2 += MinGain;
                for (N = MinN1; N != MinN2; N = N->OldSuc) {
                    if (MinForward) {
                        (N->Next = N->OldSuc)->Prev = N;
                        N->NextCost = N->Cost;
                    } else {
                        (N->Prev = N->OldSuc)->Next = N;
                        N->Prev->NextCost = N->Cost;
                    }
                }
                if (MinForward)
                    MinN2->Prev = N->OldPred;
                else {
                    MinN2->Next = N->OldPred;
                    MinN2->NextCost = N->OldPred->Cost;
                }
            }
            First = MinForward ? MinN2 : MinN1;
        }
    } while (MinN1);

    if (Cost1 < Cost2 ? !Improved1 : Cost2 < Cost1 ? !Improved2 :
        !Improved1 || !Improved2)
        return OldCost1 / Precision;

    /* Expand the best tour into a full tour */
    N = FirstNode;
    do
        N->Mark = 0;
    while ((N = N->Suc) != FirstNode);
    N = First;
    N->Mark = N;
    do {
        if (!N->Suc->Mark && (!N->V || !N->Suc->V))
            N->OldSuc = N->Suc;
        else if (!N->Pred->Mark && (!N->V || !N->Pred->V))
            N->OldSuc = N->Pred;
        else if (Cost1 <= Cost2) {
            if (N->OldSuc->Mark)
                N->OldSuc = !N->OldPred->Mark ? N->OldPred : First;
        } else if (!N->Next->Mark)
            N->OldSuc = N->Next;
        else if (!N->Prev->Mark)
            N->OldSuc = N->Prev;
        else
            N->OldSuc = First;
        N->Mark = N;
    } while ((N = N->OldSuc) != First);
    Hash = 0;
    do {
        N->OldSuc->Pred = N;
        Hash ^= Rand[N->Id] * Rand[N->OldSuc->Id];
    }
    while ((N = N->Suc = N->OldSuc) != First);
    if (TraceLevel >= 2)
        printff("IPT: " GainFormat "\n",
                (Cost1 <= Cost2 ? Cost1 : Cost2) / Precision);
    return (Cost1 <= Cost2 ? Cost1 : Cost2) / Precision;
}
//################################### MergeWithTourIPT.c end ###################################

//################################### Minimum1TreeCost.c begin ###################################
#include "LKH.h"

/*
 * The Minimum1TreeCost function returns the cost of a minimum 1-tree.
 *
 * The minimum 1-tre is found by determining the minimum spanning tree and 
 * then adding an edge corresponding to the second nearest neighbor of one 
 * of the leaves of the tree (any node which has degree 1). The leaf chosen
 * is the one that has the longest second nearest neighbor distance.
 *
 * The V-value of a node is its degree minus 2. Therefore, Norm being the 
 * sum of squares of all V-values, is a measure of a minimum 1-tree/s 
 * discrepancy from a tour. If Norm is zero, then the 1-tree constitutes a 
 * tour, and an optimal tour has been found.
 */

GainType Minimum1TreeCost(int Sparse)
{
    Node *N, *N1 = 0;
    GainType Sum = 0;
    int Max = INT_MIN;

    MinimumSpanningTree(Sparse);
    N = FirstNode;
    do {
        N->V = -2;
        Sum += N->Pi;
    }
    while ((N = N->Suc) != FirstNode);
    Sum *= -2;
    while ((N = N->Suc) != FirstNode) {
        N->V++;
        N->Dad->V++;
        Sum += N->Cost;
        N->Next = 0;
    }
    FirstNode->Dad = FirstNode->Suc;
    FirstNode->Cost = FirstNode->Suc->Cost;
    do {
        if (N->V == -1) {
            Connect(N, Max, Sparse);
            if (N->NextCost > Max && N->Next) {
                N1 = N;
                Max = N->NextCost;
            }
        }
    }
    while ((N = N->Suc) != FirstNode);
    assert(N1);
    N1->Next->V++;
    N1->V++;
    Sum += N1->NextCost;
    Norm = 0;
    do
        Norm += N->V * N->V;
    while ((N = N->Suc) != FirstNode);
    if (N1 == FirstNode)
        N1->Suc->Dad = 0;
    else {
        FirstNode->Dad = 0;
        Precede(N1, FirstNode);
        FirstNode = N1;
    }
    if (Norm == 0) {
        for (N = FirstNode->Dad; N; N1 = N, N = N->Dad)
            Follow(N, N1);
        for (N = FirstNode->Suc; N != FirstNode; N = N->Suc) {
            N->Dad = N->Pred;
            N->Cost = D(N, N->Dad);
        }
        FirstNode->Suc->Dad = 0;
    }
    return Sum;
}
//################################### Minimum1TreeCost.c end ###################################

//################################### MinimumSpanningTree.c begin ###################################
#include "LKH.h"
#include "Heap.h"

/*
 * The MinimumSpanningTree function determines a minimum spanning tree using 
 * Prim's algorithm.
 *
 * At return the Dad field of each node contains the father of the node, and 
 * the Cost field contains cost of the corresponding edge. The nodes are 
 * placed in a topological ordered list, i.e., for any node its father precedes 
 * the node in the list. The fields Pred and Suc of a node are pointers to the 
 * predecessor and successor node in this list.
 *
 * The function can be used to determine a minimum spanning tree in a dense 
 * graph, or in a sparse graph (a graph determined by a candidate set).
 *
 * When the graph is sparse a priority queue, implemented as a binary heap, 
 * is used  to speed up the determination of which edge to include next into 
 * the tree. The Rank field of a node is used to contain its priority (usually 
 * equal to the shortest distance (Cost) to nodes of the tree).        
 */

void MinimumSpanningTree(int Sparse)
{
    Node *Blue;         /* Points to the last node included in the tree */
    Node *NextBlue = 0; /* Points to the provisional next node to be included */
    Node *N;
    Candidate *NBlue;
    int d;

    Blue = N = FirstNode;
    Blue->Dad = 0;      /* The root of the tree has no father */
    if (Sparse && Blue->CandidateSet) {
        /* The graph is sparse */
        /* Insert all nodes in the heap */
        Blue->Loc = 0;  /* A blue node is not in the heap */
        while ((N = N->Suc) != FirstNode) {
            N->Dad = Blue;
            N->Cost = N->Rank = INT_MAX;
            HeapLazyInsert(N);
        }
        /* Update all neighbors to the blue node */
        for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
            if (FixedOrCommon(Blue, N)) {
                N->Dad = Blue;
                N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                N->Rank = INT_MIN;
                HeapSiftUp(N);
            } else if (!Blue->FixedTo2 && !N->FixedTo2) {
                N->Dad = Blue;
                N->Cost = N->Rank = NBlue->Cost + Blue->Pi + N->Pi;
                HeapSiftUp(N);
            }
        }
        /* Loop as long as there are more nodes to include in the tree */
        while ((NextBlue = HeapDeleteMin())) {
            Follow(NextBlue, Blue);
            Blue = NextBlue;
            /* Update all neighbors to the blue node */
            for (NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++) {
                if (!N->Loc)
                    continue;
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = NBlue->Cost + Blue->Pi + N->Pi;
                    N->Rank = INT_MIN;
                    HeapSiftUp(N);
                } else if (!Blue->FixedTo2 && !N->FixedTo2 &&
                           (d =
                            NBlue->Cost + Blue->Pi + N->Pi) < N->Cost) {
                    N->Dad = Blue;
                    N->Cost = N->Rank = d;
                    HeapSiftUp(N);
                }
            }
        }
    } else {
        /* The graph is dense */
        while ((N = N->Suc) != FirstNode)
            N->Cost = INT_MAX;
        /* Loop as long as there a more nodes to include in the tree */
        while ((N = Blue->Suc) != FirstNode) {
            int Min = INT_MAX;
            /* Update all non-blue nodes (the successors of Blue in the list) */
            do {
                if (FixedOrCommon(Blue, N)) {
                    N->Dad = Blue;
                    N->Cost = D(Blue, N);
                    NextBlue = N;
                    Min = INT_MIN;
                } else {
                    if (!Blue->FixedTo2 && !N->FixedTo2 &&
                        !Forbidden(Blue, N) &&
                        (!c || c(Blue, N) < N->Cost) &&
                        (d = D(Blue, N)) < N->Cost) {
                        N->Cost = d;
                        N->Dad = Blue;
                    }
                    if (N->Cost < Min) {
                        Min = N->Cost;
                        NextBlue = N;
                    }
                }
            }
            while ((N = N->Suc) != FirstNode);
            Follow(NextBlue, Blue);
            Blue = NextBlue;
        }
    }
}
//################################### MinimumSpanningTree.c end ###################################

//################################### NormalizeNodeList.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The NormalizeNodeList function is used to swap the Suc and Pred fields 
 * of nodes in such a way that the list of nodes constitutes a cyclic 
 * two-way list. 
 *
 * A call of the function corrupts the segment list representation.   
 */

void NormalizeNodeList()
{
    Node *t1, *t2;

    t1 = FirstNode;
    do {
        t2 = SUC(t1);
        t1->Pred = PRED(t1);
        t1->Suc = t2;
        t1->Parent = 0;
    }
    while ((t1 = t2) != FirstNode);
    Reversed = 0;
}
//################################### NormalizeNodeList.c end ###################################

//################################### NormalizeSegmentList.c begin ###################################
#include "LKH.h"

/*
 * The NormalizeSegmentList function is used to swap the Suc and Pred fields 
 * of segments in such a way that the list of segments constitutes a cyclic 
 * two-way list. 
 *
 * A call of the function corrupts the tree representation of the tour.   
 *
 * The function is called from LinKernighan.   
 */

void NormalizeSegmentList()
{
    Segment *s1, *s2;

    s1 = FirstSegment;
    do {
        if (!s1->Parent->Reversed)
            s2 = s1->Suc;
        else {
            s2 = s1->Pred;
            s1->Pred = s1->Suc;
            s1->Suc = s2;
        }
    }
    while ((s1 = s2) != FirstSegment);
}
//################################### NormalizeSegmentList.c end ###################################

//################################### OrderCandidateSet.c begin ###################################
#include "LKH.h"

/*
 * The OrderCandidateSet function augments the candidate set by using 
 * transitive relatations in the following way. If the edges (i,j) and (j,k) 
 * are contained the candidate set, then the edge (i,k) is added to the 
 * candidate set. The alpha-value of each candidate edge is computed, and 
 * the candidate edges associated with each node are ordered according to 
 * their Alpha-values.
 *   
 * The parameter MaxCandidates specifies the maximum number of candidate 
 * edges allowed for each node, and MaxAlpha puts an upper limit on their 
 * Alpha-values.
 *
 * A non-zero value of Symmetric specifies that the candidate set is to be
 * complemented such that every candidate edge is associated with both its 
 * two end nodes (in this way MaxCandidates may be exceeded).
 */

#define Ancestor OldPred        /* Nearest possible least ancestor */
#define AncestorSon OldSuc      /* Nearest son of Ancestor         */
#define OriginalCandidates Sons /* Number of original candidates   */
#define AlphaComputed OldSucExcluded    /* Have Alpha-values been computed? */
#define Level V

static int BetaValue(Node * From, Node * To);
static Candidate *FindCandidate(Node * From, Node * To);
#undef max
static int max(const int a, const int b);

void OrderCandidateSet(int MaxCandidates, GainType MaxAlpha, int Symmetric)
{
    Node *From, *To, *N;
    Candidate *NFrom, *NN;
    int Alpha, Beta;

    if (TraceLevel >= 2)
        printff("Ordering candidates ... ");
    if (MaxAlpha < 0 || MaxAlpha > INT_MAX)
        MaxAlpha = INT_MAX;
    /* Add edges from the 1-tree to the candidate set */
    if (MaxCandidates > 0) {
        From = FirstNode;
        do {
            if ((To = From->Dad)) {
                AddCandidate(From, To, From->Cost, 0);
                AddCandidate(To, From, From->Cost, 0);
            }
        }
        while ((From = From->Suc) != FirstNode);
        AddCandidate(FirstNode, FirstNode->Next, FirstNode->NextCost, 0);
        AddCandidate(FirstNode->Next, FirstNode, FirstNode->NextCost, 0);
    }

    From = FirstNode;
    do {
        From->AlphaComputed = 0;
        From->Sons = 0;
    } while ((From = From->Suc) != FirstNode);
    From = FirstNode->Suc;
    From->Level = 0;
    From->Ancestor = From->AncestorSon = From;
    From->Beta = INT_MIN;
    From = From->Suc;
    do {
        From->Ancestor = To = From->Dad;
        From->Beta = !FixedOrCommon(From, To) ? From->Cost : INT_MIN;
        From->Level = To->Level + 1;
        To->Sons++;
        From->AncestorSon = From;
    }
    while ((From = From->Suc) != FirstNode);

    From = FirstNode->Suc->Suc;
    do {
        To = From->Dad;
        if (To->Sons == 1) {
            From->Beta = max(To->Beta, From->Beta);
            From->Ancestor = To->Ancestor;
            From->AncestorSon = To->AncestorSon;
        }
    }
    while ((From = From->Suc) != FirstNode);

    /* Compute Alpha-values for candidates */
    do {
        for (NFrom = From->CandidateSet; NFrom && (To = NFrom->To);
             NFrom++) {
            if (FixedOrCommon(From, To))
                NFrom->Alpha = INT_MIN;
            else if (From->FixedTo2 || To->FixedTo2 || Forbidden(From, To))
                NFrom->Alpha = INT_MAX;
            else if (To->AlphaComputed && (NN = FindCandidate(To, From)))
                NFrom->Alpha = NN->Alpha;
            else {
                Beta = BetaValue(From, To);
                NFrom->Alpha =
                    Beta != INT_MIN ? max(NFrom->Cost - Beta, 0) : INT_MAX;
                if (NFrom->Alpha > MaxAlpha &&
                    !DelaunayPure && CandidateSetType == DELAUNAY)
                    NFrom->Alpha = INT_MAX;
            }
        }
        From->AlphaComputed = 1;
    }
    while ((From = From->Suc) != FirstNode);

    if (MaxCandidates > 0 && !DelaunayPure && CandidateSetType == DELAUNAY) {
        do {
            int Count = 0;
            for (NFrom = From->CandidateSet; NFrom->To; NFrom++)
                Count++;
            From->OriginalCandidates = Count;
            From->Mark = 0;
            From->AlphaComputed = 0;
        }
        while ((From = From->Suc) != FirstNode);

        /* Augment the original candidate set */
        do {
            int i, j;
            From->Mark = From;
            for (i = 0; i < From->OriginalCandidates; i++) {
                N = From->CandidateSet[i].To;
                N->Mark = From;
                for (j = 0; j < N->OriginalCandidates; j++) {
                    To = N->CandidateSet[j].To;
                    if (To->Mark == From)
                        continue;
                    To->Mark = From;
                    if (FindCandidate(From, To))
                        continue;
                    if (FixedOrCommon(From, To))
                        Alpha = INT_MIN;
                    else if (From->FixedTo2 || To->FixedTo2 ||
                             Forbidden(From, To))
                        continue;
                    else if ((NN = FindCandidate(To, From)))
                        Alpha = NN->Alpha;
                    else {
                        if (To->AlphaComputed)
                            continue;
                        Beta = BetaValue(From, To);
                        if (Beta == INT_MIN)
                            continue;
                        Alpha = max(D(From, To) - Beta, 0);
                    }
                    if (Alpha <= MaxAlpha)
                        AddCandidate(From, To, D(From, To), Alpha);
                }
            }
            From->AlphaComputed = 1;
        }
        while ((From = From->Suc) != FirstNode);
    }

    /* Order candidates according to their Alpha-values */
    ResetCandidateSet();
    if (MaxCandidates > 0)
        TrimCandidateSet(MaxCandidates);
    AddTourCandidates();
    if (Symmetric)
        SymmetrizeCandidateSet();
    if (TraceLevel >= 2)
        printff("done\n");
}

/*
 * The BetaValue function computes the largest edge cost on the path 
 * between two given nodes in the minimum spanning tree.    
 */

static int BetaValue(Node * From, Node * To)
{
    Node *N1 = From, *N2 = To;
    int Beta = INT_MIN;

    if (To == From->Dad)
        return From->Cost;
    if (From == To->Dad)
        return To->Cost;
    if (From == FirstNode || To == FirstNode)
        return FirstNode->NextCost;

    /* Go upwards in the tree until the least common ancestor is met */
    while (N1->Ancestor != N2->Ancestor) {
        if (N1->Level > N2->Level) {
            if (N1->Beta > Beta)
                Beta = N1->Beta;
            N1 = N1->Ancestor;
        } else {
            if (N2->Beta > Beta)
                Beta = N2->Beta;
            N2 = N2->Ancestor;
        }
    }
    if (N1 == N2)
        return Beta;
    if (N1->AncestorSon != N2->AncestorSon)
        return max(Beta, max(N1->Beta, N2->Beta));
    if (N1->Level < N2->Level) {
        Node *t = N1;
        N1 = N2;
        N2 = t;
    }
    if (N1->Beta > N2->Beta)
        return max(Beta, N1->Beta);
    while (N1 != N2) {
        if (N1->Cost > Beta)
            Beta = N1->Cost;
        N1 = N1->Dad;
    }
    return Beta;
}

/*
 * The FindCandidate function returns the Candidate structure that is
 * associated with the node From and is pointing to the node To. The
 * function returns 0 if the search fails.     
 */

static Candidate *FindCandidate(Node * From, Node * To)
{
    Candidate *NFrom;
    for (NFrom = From->CandidateSet; NFrom->To; NFrom++)
        if (NFrom->To == To)
            return NFrom;
    return 0;
}

static int max(const int a, const int b)
{
    return a > b ? a : b;
}
//################################### OrderCandidateSet.c end ###################################

//################################### PatchCycles.c begin ###################################
#include "Segment.h"
#include "LKH.h"
#include "Sequence.h"

/*
 * The PatchCycles function attempts to improve the tour by patching the
 * M >= 2 cycles that would appear if the move defined by t[1..2k] and
 * incl[1..2k] was made. If the composite move results in a shorter
 * tour, then the move is made, and the function returns the gain.
 *
 * On entry, Gain is the gain that could be obtained by making the non-
 * feasible move defined by t[1..2k] and incl[1..2k].
 *
 * The function tries to patch the cycles by interleaving the alternating
 * path represented by t with one or more alternating cycles.
 *
 * The function is called from BestKOptMove.
 */

static GainType PatchCyclesRec(int k, int m, int M, GainType G0);
static int ShortestCycle(int M, int k);
static int Cycle(Node * N, int k);

static int CurrentCycle, Patchwork = 0, RecLevel = 0;
#define MaxPatchwork Dimension

/*
 * The PatchCycles function tries to find a gainful move by patching the
 * cycles that would occur if the move represented by t[1..2k] and incl[1..2k]
 * was made using one one or more alternating cycles.
 * The alternating cycles are put in continuation of t, starting at 2k+1.
 */

GainType PatchCycles(int k, GainType Gain)
{
    Node *s1, *s2, *sStart, *sStop;
    GainType NewGain;
    int M, i;

    FindPermutation(k);
    M = Cycles(k);
    if (M == 1 && Gain > 0) {
        MakeKOptMove(k);
        return Gain;
    }
    if (M == 1 || M > PatchingC || k + M > NonsequentialMoveType)
        return 0;
    if (RecLevel == 0)
        Patchwork = 0;
    CurrentCycle = ShortestCycle(M, k);
    for (i = 0; i < k; i++) {
        if (cycle[p[2 * i]] != CurrentCycle)
            continue;
        sStart = t[p[2 * i]];
        sStop = t[p[2 * i + 1]];
        for (s1 = sStart; s1 != sStop; s1 = s2) {
            s2 = SUC(s1);
            if (FixedOrCommon(s1, s2))
                continue;
            if (++Patchwork > MaxPatchwork)
                return 0;
            t[2 * k + 1] = s1;
            t[2 * k + 2] = s2;
            MarkDeleted(s1, s2);
            /* Find a set of gainful alternating cycles */
            NewGain = PatchCyclesRec(k, 2, M, Gain + C(s1, s2));
            UnmarkDeleted(s1, s2);
            if (NewGain > 0)
                return NewGain;
        }
    }
    return 0;
}

static GainType PatchCyclesRec(int k, int m, int M, GainType G0)
{
    Node *s1, *s2, *s3, *s4, *s5, *s6, *S3 = 0, *S4 = 0;
    Candidate *Ns2, *Ns4;
    GainType G1, G2, G3, G4, Gain, CloseUpGain,
        BestCloseUpGain = PatchingAExtended ? MINUS_INFINITY : 0;
    int X4, X6;
    int i, NewCycle, *cycleSaved = 0, *pSaved = 0;
    int Breadth2 = 0, Breadth4;

    s1 = t[2 * k + 1];
    s2 = t[i = 2 * (k + m) - 2];
    incl[incl[i] = i + 1] = i;

    /* Choose (s2,s3) as a candidate edge emanating from s2 */
    for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
        if (s3 == s2->Pred || s3 == s2->Suc || Added(s2, s3) ||
            (NewCycle = Cycle(s3, k)) == CurrentCycle)
            continue;
        if (++Breadth2 > MaxBreadth)
            break;
        MarkAdded(s2, s3);
        t[2 * (k + m) - 1] = s3;
        G1 = G0 - Ns2->Cost;
        /* Choose s4 as one of s3's two neighbors on the tour */
        for (X4 = 1; X4 <= 2; X4++) {
            s4 = X4 == 1 ? s3->Pred : s3->Suc;
            if (FixedOrCommon(s3, s4) || Deleted(s3, s4))
                continue;
            MarkDeleted(s3, s4);
            t[2 * (k + m)] = s4;
            G2 = G1 + C(s3, s4);
            if (M > 2) {
                if (!cycleSaved) {
                    cycleSaved = (int *) malloc(2 * k * sizeof(int));
                    memcpy(cycleSaved, cycle + 1, 2 * k * sizeof(int));
                }
                for (i = 1; i <= 2 * k; i++)
                    if (cycle[i] == NewCycle)
                        cycle[i] = CurrentCycle;
                /* Extend the current alternating path */
                if ((Gain = PatchCyclesRec(k, m + 1, M - 1, G2)) > 0) {
                    UnmarkAdded(s2, s3);
                    UnmarkDeleted(s3, s4);
                    goto End_PatchCyclesRec;
                }
                memcpy(cycle + 1, cycleSaved, 2 * k * sizeof(int));
                if (PatchingA >= 2 && Patchwork < MaxPatchwork &&
                    k + M < NonsequentialMoveType &&
                    !Forbidden(s4, s1) &&
                    (!PatchingARestricted || IsCandidate(s4, s1))) {
                    GainType Bound = BestCloseUpGain >= 0 ||
                        IsCandidate(s4, s1) ? BestCloseUpGain : 0;
                    if ((!c || G2 - c(s4, s1) > Bound) &&
                        (CloseUpGain = G2 - C(s4, s1)) > Bound) {
                        S3 = s3;
                        S4 = s4;
                        BestCloseUpGain = CloseUpGain;
                    }
                }
            } else if (!Forbidden(s4, s1) && (!c || G2 - c(s4, s1) > 0)
                       && (Gain = G2 - C(s4, s1)) > 0) {
                incl[incl[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
                MakeKOptMove(k + m);
                UnmarkAdded(s2, s3);
                UnmarkDeleted(s3, s4);
                goto End_PatchCyclesRec;
            }
            UnmarkDeleted(s3, s4);
        }
        UnmarkAdded(s2, s3);
    }
    if (M == 2 && !PatchingCRestricted) {
        /* Try to patch the two cycles by a sequential 3-opt move */
        incl[incl[2 * (k + m)] = 2 * (k + m) + 1] = 2 * (k + m);
        incl[incl[2 * k + 1] = 2 * (k + m) + 2] = 2 * k + 1;
        Breadth2 = 0;
        /* Choose (s2,s3) as a candidate edge emanating from s2 */
        for (Ns2 = s2->CandidateSet; (s3 = Ns2->To); Ns2++) {
            if (s3 == s2->Pred || s3 == s2->Suc || Added(s2, s3))
                continue;
            if (++Breadth2 > MaxBreadth)
                break;
            t[2 * (k + m) - 1] = s3;
            G1 = G0 - Ns2->Cost;
            NewCycle = Cycle(s3, k);
            /* Choose s4 as one of s3's two neighbors on the tour */
            for (X4 = 1; X4 <= 2; X4++) {
                s4 = X4 == 1 ? s3->Pred : s3->Suc;
                if (FixedOrCommon(s3, s4) || Deleted(s3, s4))
                    continue;
                t[2 * (k + m)] = s4;
                G2 = G1 + C(s3, s4);
                Breadth4 = 0;
                /* Choose (s4,s5) as a candidate edge emanating from s4 */
                for (Ns4 = s4->CandidateSet; (s5 = Ns4->To); Ns4++) {
                    if (s5 == s4->Pred || s5 == s4->Suc || s5 == s1 ||
                        Added(s4, s5) ||
                        (NewCycle == CurrentCycle &&
                         Cycle(s5, k) == CurrentCycle))
                        continue;
                    if (++Breadth4 > MaxBreadth)
                        break;
                    G3 = G2 - Ns4->Cost;
                    /* Choose s6 as one of s5's two neighbors on the tour */
                    for (X6 = 1; X6 <= 2; X6++) {
                        s6 = X6 == 1 ? s5->Pred : s5->Suc;
                        if (s6 == s1 || Forbidden(s6, s1)
                            || FixedOrCommon(s5, s6)
                            || Deleted(s5, s6)
                            || Added(s6, s1))
                            continue;
                        G4 = G3 + C(s5, s6);
                        if ((!c || G4 - c(s6, s1) > 0) &&
                            (Gain = G4 - C(s6, s1)) > 0) {
                            if (!pSaved) {
                                pSaved = (int *) malloc(2 * k * sizeof(int));
                                memcpy(pSaved, p + 1, 2 * k * sizeof(int));
                            }
                            t[2 * (k + m) + 1] = s5;
                            t[2 * (k + m) + 2] = s6;
                            if (FeasibleKOptMove(k + m + 1)) {
                                MakeKOptMove(k + m + 1);
                                goto End_PatchCyclesRec;
                            }
                            memcpy(p + 1, pSaved, 2 * k * sizeof(int));
                            for (i = 1; i <= 2 * k; i++)
                                q[p[i]] = i;
                        }
                    }
                }
            }
        }
    }
    Gain = 0;
    if (S4) {
        int OldCycle = CurrentCycle;
        if (!pSaved) {
            pSaved = (int *) malloc(2 * k * sizeof(int));
            memcpy(pSaved, p + 1, 2 * k * sizeof(int));
        }
        t[2 * (k + m) - 1] = S3;
        t[2 * (k + m)] = S4;
        incl[incl[2 * k + 1] = 2 * (k + m)] = 2 * k + 1;
        /* Find a new alternating cycle */
        PatchingA--;
        RecLevel++;
        MarkAdded(s2, S3);
        MarkDeleted(S3, S4);
        MarkAdded(S4, s1);
        Gain = PatchCycles(k + m, BestCloseUpGain);
        UnmarkAdded(s2, S3);
        UnmarkDeleted(S3, S4);
        UnmarkAdded(S4, s1);
        RecLevel--;
        PatchingA++;
        if (Gain <= 0) {
            memcpy(cycle + 1, cycleSaved, 2 * k * sizeof(int));
            memcpy(p + 1, pSaved, 2 * k * sizeof(int));
            for (i = 1; i <= 2 * k; i++)
                q[p[i]] = i;
            CurrentCycle = OldCycle;
        }
    }

  End_PatchCyclesRec:
    free(cycleSaved);
    free(pSaved);
    return Gain;
}

/*
 * The Cycle function returns the number of the cycle containing
 * a given node, N.
 *
 * Time complexity: O(log k).
 */

static int Cycle(Node * N, int k)
{
    /* Binary search */
    int Low = 1, High = k;
    while (Low < High) {
        int Mid = (Low + High) / 2;
        if (BETWEEN(t[p[2 * Low]], N, t[p[2 * Mid + 1]]))
            High = Mid;
        else
            Low = Mid + 1;
    }
    return cycle[p[2 * Low]];
}

/*
 * The ShortestCycle function returns the number of the cycle with
 * the smallest number of nodes. Note however that if the two-level
 * list is used, the number of nodes of each cycle is only approximate
 * (for efficiency reasons).
 *
 * Time complexity: O(k + M), where M = Cycles(k).
 *
 * The function may only be called after a call of the Cycles function.
 */

static int ShortestCycle(int M, int k)
{
    int i, Cycle, MinCycle = 0;
    int *Size, MinSize = INT_MAX;

    Size = (int *) calloc(1 + M, sizeof(int));
    p[0] = p[2 * k];
    for (i = 0; i < 2 * k; i += 2)
        Size[cycle[p[i]]] += SegmentSize(t[p[i]], t[p[i + 1]]);
    for (Cycle = 1; Cycle <= M; Cycle++) {
        if (Size[Cycle] < MinSize) {
            MinSize = Size[Cycle];
            MinCycle = Cycle;
        }
    }
    free(Size);
    return MinCycle;
}
//################################### PatchCycles.c end ###################################

//################################### printff.c begin ###################################
#include <stdio.h>
#include <stdarg.h>

/* 
 * The printff function prints a message and flushes stdout.
 */

void printff(const char *fmt, ...)
{
    va_list args;

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
}
//################################### printff.c end ###################################

//################################### PrintParameters.c begin ###################################
#include "LKH.h"
#include "Genetic.h"

/*
 * The PrintParameters function prints the problem parameters to 
 * standard output. 
*/

void PrintParameters()
{
    int i;

    printff("ASCENT_CANDIDATES = %d\n", AscentCandidates);
    printff("BACKBONE_TRIALS = %d\n", BackboneTrials);
    printff("BACKTRACKING = %s\n", Backtracking ? "YES" : "NO");
    if (CandidateFiles == 0)
        printff("# CANDIDATE_FILE =\n");
    else
        for (i = 0; i < CandidateFiles; i++)
            printff("CANDIDATE_FILE = %s\n", CandidateFileName[i]);
    printff("CANDIDATE_SET_TYPE = %s%s\n",
            CandidateSetType == ALPHA ? "ALPHA" :
            CandidateSetType == DELAUNAY ? "DELAUNAY" :
            CandidateSetType == NN ? "NEAREST-NEIGHBOR" :
            CandidateSetType == POPMUSIC ? "POPMUSIC" :
            CandidateSetType == QUADRANT ? "QUADRANT" : "",
            DelaunayPure ? " PURE" : "");
    if (EdgeFiles == 0)
        printff("# EDGE_FILE =\n");
    else
        for (i = 0; i < EdgeFiles; i++)
            printff("EDGE_FILE = %s\n", EdgeFileName[i]);
    if (Excess >= 0)
        printff("EXCESS = %g\n", Excess);
    else
        printff("# EXCESS =\n");
    printff("EXTRA_CANDIDATES = %d %s\n",
            ExtraCandidates,
            ExtraCandidateSetSymmetric ? "SYMMETRIC" : "");
    printff("EXTRA_CANDIDATE_SET_TYPE = %s\n",
            ExtraCandidateSetType == NN ? "NEAREST-NEIGHBOR" :
            ExtraCandidateSetType == POPMUSIC ? "POPMUSIC" :
            ExtraCandidateSetType == QUADRANT ? "QUADRANT" : "");
    printff("GAIN23 = %s\n", Gain23Used ? "YES" : "NO");
    printff("GAIN_CRITERION = %s\n", GainCriterionUsed ? "YES" : "NO");
    if (InitialPeriod >= 0)
        printff("INITIAL_PERIOD = %d\n", InitialPeriod);
    else
        printff("# INITIAL_PERIOD =\n");
    printff("INITIAL_STEP_SIZE = %d\n", InitialStepSize);
    printff("INITIAL_TOUR_ALGORITHM = %s\n",
            InitialTourAlgorithm == BORUVKA ? "BORUVKA" :
            InitialTourAlgorithm == GREEDY ? "GREEDY" :
            InitialTourAlgorithm == MOORE ? "MOORE" :
            InitialTourAlgorithm == NEAREST_NEIGHBOR ? "NEAREST-NEIGHBOR" :
            InitialTourAlgorithm == QUICK_BORUVKA ? "QUICK-BORUVKA" :
            InitialTourAlgorithm == SIERPINSKI ? "SIERPINSKI" : "WALK");
    printff("%sINITIAL_TOUR_FILE = %s\n",
            InitialTourFileName ? "" : "# ",
            InitialTourFileName ? InitialTourFileName : "");
    printff("INITIAL_TOUR_FRACTION = %0.3f\n", InitialTourFraction);
    printff("%sINPUT_TOUR_FILE = %s\n",
            InputTourFileName ? "" : "# ",
            InputTourFileName ? InputTourFileName : "");
    printff("KICK_TYPE = %d\n", KickType);
    printff("KICKS = %d\n", Kicks);
    if (MaxBreadth == INT_MAX)
        printff("# MAX_BREADTH =\n");
    else
        printff("MAX_BREADTH = %d\n", MaxBreadth);
    printff("MAX_CANDIDATES = %d %s\n",
            MaxCandidates, CandidateSetSymmetric ? "SYMMETRIC" : "");
    if (MaxSwaps >= 0)
        printff("MAX_SWAPS = %d\n", MaxSwaps);
    else
        printff("# MAX_SWAPS =\n");
    if (MaxTrials >= 0)
        printff("MAX_TRIALS = %d\n", MaxTrials);
    else
        printff("# MAX_TRIALS =\n");
    if (MergeTourFiles == 0)
        printff("# MERGE_TOUR_FILE =\n");
    else
        for (i = 0; i < MergeTourFiles; i++)
            printff("MERGE_TOUR_FILE = %s\n", MergeTourFileName[i]);
    printff("MOVE_TYPE = %d\n", MoveType);
    printff("%sNONSEQUENTIAL_MOVE_TYPE = %d\n",
            PatchingA > 1 ? "" : "# ", NonsequentialMoveType);
    if (Optimum == MINUS_INFINITY)
        printff("# OPTIMUM =\n");
    else
        printff("OPTIMUM = " GainFormat "\n", Optimum);
    printff("%sOUTPUT_TOUR_FILE = %s\n",
            OutputTourFileName ? "" : "# ",
            OutputTourFileName ? OutputTourFileName : "");
    printff("PATCHING_A = %d %s\n", PatchingA,
            PatchingARestricted ? "RESTRICTED" :
            PatchingAExtended ? "EXTENDED" : "");
    printff("PATCHING_C = %d %s\n", PatchingC,
            PatchingCRestricted ? "RESTRICTED" :
            PatchingCExtended ? "EXTENDED" : "");
    printff("%sPI_FILE = %s\n",
            PiFileName ? "" : "# ", PiFileName ? PiFileName : "");
    printff("POPMUSIC_INITIAL_TOUR = %s\n",
            POPMUSIC_InitialTour ? "YES" : "NO");
    printff("POPMUSIC_MAX_NEIGHBORS = %d\n", POPMUSIC_MaxNeighbors);
    printff("POPMUSIC_SAMPLE_SIZE = %d\n", POPMUSIC_SampleSize);
    printff("POPMUSIC_SOLUTIONS = %d\n", POPMUSIC_Solutions);
    printff("POPMUSIC_TRIALS = %d\n", POPMUSIC_Trials);
    if (MaxPopulationSize == 0)
        printff("# ");
    printff("POPULATION_SIZE = %d\n", MaxPopulationSize);
    printff("PRECISION = %d\n", Precision);
    printff("%sPROBLEM_FILE = %s\n",
            ProblemFileName ? "" : "# ",
            ProblemFileName ? ProblemFileName : "");
   printff("RECOMBINATION = %s\n",
           Recombination == IPT ? "IPT" :
           Recombination == GPX2 ? "GPX2" :
           Recombination == CLARIST ? "CLARIST" :
           "UNKNOWN");
    printff("RESTRICTED_SEARCH = %s\n", RestrictedSearch ? "YES" : "NO");
    printff("RUNS = %d\n", Runs);
    printff("SEED = %u\n", Seed);
    printff("STOP_AT_OPTIMUM = %s\n", StopAtOptimum ? "YES" : "NO");
    printff("SUBGRADIENT = %s\n", Subgradient ? "YES" : "NO");
    if (SubproblemSize == 0)
        printff("# SUBPROBLEM_SIZE =\n");
    else
        printff("SUBPROBLEM_SIZE = %d%s%s%s\n", SubproblemSize,
                DelaunayPartitioning ? " DELAUNAY" :
                KarpPartitioning ? " KARP" :
                KCenterPartitioning ? " K-CENTER" :
                KMeansPartitioning ? " K-MEANS" :
                MoorePartitioning ? " MOORE" :
                RohePartitioning ? " ROHE" :
                SierpinskiPartitioning ? " SIERPINSKI" : "",
                SubproblemBorders ? " BORDERS" : "",
                SubproblemsCompressed ? " COMPRESSED" : "");
    printff("%sSUBPROBLEM_TOUR_FILE = %s\n",
            SubproblemTourFileName ? "" : "# ",
            SubproblemTourFileName ? SubproblemTourFileName : "");
    printff("SUBSEQUENT_MOVE_TYPE = %d\n",
            SubsequentMoveType == 0 ? MoveType : SubsequentMoveType);
    printff("SUBSEQUENT_PATCHING = %s\n",
            SubsequentPatching ? "YES" : "NO");
    if (TimeLimit == DBL_MAX)
        printff("# TIME_LIMIT =\n");
    else
        printff("TIME_LIMIT = %0.1f\n", TimeLimit);
    if (TotalTimeLimit == DBL_MAX)
        printff("# TOTAL_TIME_LIMIT =\n");
    else
        printff("TOTAL_TIME_LIMIT = %0.1f\n", TotalTimeLimit);
    printff("%sTOUR_FILE = %s\n",
            TourFileName ? "" : "# ", TourFileName ? TourFileName : "");
    printff("TRACE_LEVEL = %d\n\n", TraceLevel);
}
//################################### PrintParameters.c end ###################################

//################################### Random.c begin ###################################
/*
 * This file contains a portable random generator. It will give
 * identical sequences of random integers for any platform with
 * at least 32-bit integers.
 *
 * A version of this generator is described in J. Bentley's column, 
 * "The Software Exploratorium", Unix Review 1991. It is based on 
 * Algorithm A in D. E. Knuth, The Art of Computer Programming, 
 * Vol 2, Section 3.2.2, pp. 172.  
 *  
 * The Random function returns a pseudo-random integer in the range
 * 0...INT_MAX-1.
 *   
 * The SRandom function uses the given seed for a new sequence of
 * pseudo-random numbers.  
 */

unsigned Random(void);
void SRandom(unsigned Seed);

#undef STDLIB_RANDOM
/* #define STDLIB_RANDOM */

#ifdef STDLIB_RANDOM
#include <stdlib.h>
unsigned Random()
{
    return rand();
}

void SRandom(unsigned Seed)
{
    srand(Seed);
}

#else

#include <limits.h>
#define PRANDMAX INT_MAX

static int a = 0, b = 24, arr[55], initialized = 0;

unsigned Random()
{
    int t;

    if (!initialized)
        SRandom(7913);
    if (a-- == 0)
        a = 54;
    if (b-- == 0)
        b = 54;
    if ((t = arr[a] - arr[b]) < 0)
        t += PRANDMAX;
    return (arr[a] = t);
}

void SRandom(unsigned Seed)
{
    int i, ii, last, next;

    Seed %= PRANDMAX;
    arr[0] = last = Seed;
    for (next = i = 1; i < 55; i++) {
        ii = (21 * i) % 55;
        arr[ii] = next;
        if ((next = last - next) < 0)
            next += PRANDMAX;
        last = arr[ii];
    }
    initialized = 1;
    a = 0;
    b = 24;
    for (i = 0; i < 165; i++)
        Random();
}

#endif
//################################### Random.c end ###################################

//################################### RecordBestTour.c begin ###################################
#include "LKH.h"

/*
 * The RecordBestTour function records the current best tour in the BestTour 
 * array. 
 *
 * The function is called by LKHmain each time a run has resulted in a
 * shorter tour. Thus, when the predetermined number of runs have been
 * completed, BestTour contains an array representation of the best tour
 * found.    
 */

void RecordBestTour()
{
    int i;

    for (i = 0; i <= DimensionSaved; i++)
        BestTour[i] = BetterTour[i];
}
//################################### RecordBestTour.c end ###################################

//################################### RecordBetterTour.c begin ###################################
#include "LKH.h"

/*
 * The RecordBetterTour function is called by FindTour each time
 * the LinKernighan function has returned a better tour.
 *
 * The function records the tour in the BetterTour array and in the
 * BestSuc field of each node. Furthermore, for each node the previous
 * value of BestSuc is saved in the NextBestSuc field.
 *
 * Recording a better tour in the BetterTour array when the problem is
 * asymmetric requires special treatment since the number of nodes has
 * been doubled.
 */

void RecordBetterTour()
{
    Node *N = FirstNode, *Stop = N;

    // ATSP
    if (Stop->Id > DimensionSaved)
        Stop = N = Stop->Suc;
    if (N->Suc->Id != DimensionSaved + N->Id) {
        int i = 1;
        do
            if (N->Id <= DimensionSaved)
                BetterTour[i++] = N->Id;
        while ((N = N->Suc) != Stop);
    } else {
        int i = DimensionSaved;
        do
            if (N->Id <= DimensionSaved)
                BetterTour[i--] = N->Id;
        while ((N = N->Suc) != Stop);
    }
    
    BetterTour[0] = BetterTour[DimensionSaved];
    N = FirstNode;
    do {
        N->NextBestSuc = N->BestSuc;
        N->BestSuc = N->Suc;
    }
    while ((N = N->Suc) != FirstNode);
}
//################################### RecordBetterTour.c end ###################################

//################################### RemoveFirstActive.c begin ###################################
#include "LKH.h"

/* 
 * The RemoveFirstActive function removes the first node in the list 
 * of "active" nodes (i.e., nodes to be tried as an anchor node, t1,
 * by the LinKernighan algorithm).
 *
 * The function returns a pointer to the removed node. 
 *
 * The list must not be empty before the call. 
 */

Node *RemoveFirstActive()
{
    Node *N = FirstActive;
    if (FirstActive == LastActive)
        FirstActive = LastActive = 0;
    else
        LastActive->Next = FirstActive = FirstActive->Next;
    if (N)
        N->Next = 0;
    return N;
}
//################################### RemoveFirstActive.c end ###################################

//################################### ResetCandidateSet.c begin ###################################
#include "LKH.h"

/*
 * Each time a trial has resulted in a shorter tour the candidate set is
 * adjusted (by AdjustCandidateSet). The ResetCandidates function resets
 * the candidate set. The original order is re-established (using, and 
 * edges with Alpha == INT_MAX are excluded.
 *
 * The function is called from FindTour and OrderCandidates.
 */

void ResetCandidateSet()
{
    Node *From;
    Candidate *NFrom, *NN, Temp;

    From = FirstNode;
    /* Loop for all nodes */
    do {
        if (!From->CandidateSet)
            continue;
        /* Reorder the candidate array of From */
        for (NFrom = From->CandidateSet; NFrom->To; NFrom++) {
            Temp = *NFrom;
            for (NN = NFrom - 1;
                 NN >= From->CandidateSet &&
                 (Temp.Alpha < NN->Alpha ||
                  (Temp.Alpha == NN->Alpha && Temp.Cost < NN->Cost)); NN--)
                *(NN + 1) = *NN;
            *(NN + 1) = Temp;
        }
        NFrom--;
        /* Remove included edges */
        while (NFrom >= From->CandidateSet + 2 && NFrom->Alpha == INT_MAX)
            NFrom--;
        NFrom++;
        NFrom->To = 0;
        /* Remove impossible candidates */
        for (NFrom = From->CandidateSet; NFrom->To; NFrom++) {
            if (!IsPossibleCandidate(From, NFrom->To)) {
                for (NN = NFrom; NN->To; NN++)
                    *NN = *(NN + 1);
                NFrom--;
            }
        }
    }
    while ((From = From->Suc) != FirstNode);
}
//################################### ResetCandidateSet.c end ###################################

//################################### RestoreTour.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The RestoreTour function is used to undo a series of moves. The function 
 * restores the tour from SwapStack, the stack of 2-opt moves. A bad sequence 
 * of moves is undone by unstacking the 2-opt moves and making the inverse 
 * 2-opt moves in this reversed sequence.
 */

void RestoreTour()
{
    Node *t1, *t2, *t3, *t4;

    /* Loop as long as the stack is not empty */
    while (Swaps > 0) {
        /* Undo topmost 2-opt move */
        Swaps--;
        t1 = SwapStack[Swaps].t1;
        t2 = SwapStack[Swaps].t2;
        t3 = SwapStack[Swaps].t3;
        t4 = SwapStack[Swaps].t4;
        Swap1(t3, t2, t1);
        Swaps--;
        /* Make edges (t1,t2) and (t2,t3) excludable again */
        t1->OldPredExcluded = t1->OldSucExcluded = 0;
        t2->OldPredExcluded = t2->OldSucExcluded = 0;
        t3->OldPredExcluded = t3->OldSucExcluded = 0;
        t4->OldPredExcluded = t4->OldSucExcluded = 0;
    }
}
//################################### RestoreTour.c end ###################################

//################################### SegmentSize.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/*
 * The SegmentSize function returns the number of nodes in the 
 * tour segment between two given nodes in the current direction. 
 * Note, however, that if the two-level or three-level tree is used,
 * the number of nodes is only approximate (for efficiency reasons).
 * 
 * Time complexity: O(1).
 */

#ifdef ONE_LEVEL_TREE

int SegmentSize(Node * ta, Node * tb)
{
    int n = !Reversed ? tb->Rank - ta->Rank : ta->Rank - tb->Rank;
    return (n < 0 ? n + Dimension : n) + 1;
}

#elif defined TWO_LEVEL_TREE

int SegmentSize(Node * ta, Node * tb)
{
    Segment *Pa, *Pb;
    int nLeft, nMid, nRight;

    Pa = ta->Parent;
    Pb = tb->Parent;
    if (Pa == Pb) {
        int n = Reversed == Pa->Reversed ? tb->Rank - ta->Rank :
            ta->Rank - tb->Rank;
        return (n < 0 ? n + Dimension : n) + 1;
    }
    nLeft =
        Reversed ==
        Pa->Reversed ? Pa->Last->Rank - ta->Rank : ta->Rank -
        Pa->First->Rank;
    if (nLeft < 0)
        nLeft += Pa->Size;
    nMid = !Reversed ? Pb->Rank - Pa->Rank : Pa->Rank - Pb->Rank;
    if (nMid < 0)
        nMid += Groups;
    nMid = nMid == 2 ? (!Reversed ? Pa->Suc : Pa->Pred)->Size
        : (nMid - 1) * GroupSize;
    nRight =
        Reversed ==
        Pb->Reversed ? tb->Rank -
        Pb->First->Rank : Pb->Last->Rank - tb->Rank;
    if (nRight < 0)
        nRight += Pb->Size;
    return nLeft + nMid + nRight + 2;
}

#elif defined  THREE_LEVEL_TREE

int SegmentSize(Node * ta, Node * tb)
{
    Segment *Pa, *Pb;
    SSegment *PPa, *PPb;
    int n, nLeft, nMid, nRight;

    Pa = ta->Parent;
    Pb = tb->Parent;
    PPa = Pa->Parent;
    PPb = Pb->Parent;
    if (Pa == Pb) {
        n = Reversed == (Pa->Reversed != PPa->Reversed) ?
            tb->Rank - ta->Rank : ta->Rank - tb->Rank;
        if (n < 0)
            n += Dimension;
    } else if (PPa == PPb) {
        nLeft =
            Reversed == (Pa->Reversed != PPa->Reversed) ?
            Pa->Last->Rank - ta->Rank : ta->Rank - Pa->First->Rank;
        if (nLeft < 0)
            nLeft += Pa->Size;
        nMid =
            Reversed == PPa->Reversed ?
            Pb->Rank - Pa->Rank : Pa->Rank - Pb->Rank;
        if (nMid < 0)
            nMid += Groups;
        nMid = nMid == 2 ?
            (Reversed == PPa->Reversed ? Pa->Suc : Pa->Pred)->Size
            : (nMid - 1) * GroupSize;
        nRight =
            (Reversed != PPa->Reversed) ==
            (Pb->Reversed != PPb->Reversed) ? tb->Rank -
            Pb->First->Rank : Pb->Last->Rank - tb->Rank;
        if (nRight < 0)
            nRight += Pb->Size;
        n = nLeft + nMid + nRight + 1;
    } else {
        nLeft =
            Reversed == PPa->Reversed ?
            PPa->Last->Rank - Pa->Rank : Pa->Rank - PPa->First->Rank;
        if (nLeft < 0)
            nLeft += PPa->Size;
        if (nLeft > 0)
            nLeft *= GroupSize;
        else {
            nLeft =
                Reversed == (Pa->Reversed != PPa->Reversed) ?
                Pa->Last->Rank - ta->Rank : ta->Rank - Pa->First->Rank;
            if (nLeft < 0)
                nLeft += Pa->Size;
        }
        nMid = !Reversed ? PPb->Rank - PPa->Rank : PPa->Rank - PPb->Rank;
        if (nMid < 0)
            nMid += SGroups;
        nMid = nMid == 2 ?
            (!Reversed ? PPa->Suc : PPa->Pred)->Size
            : (nMid - 1) * SGroupSize;
        nMid *= GroupSize;
        nRight =
            Reversed == PPb->Reversed ? Pb->Rank -
            PPb->First->Rank : PPb->Last->Rank - Pb->Rank;
        if (nRight < 0)
            nRight += PPb->Size;
        if (nRight > 0)
            nRight *= GroupSize;
        else {
            nRight =
                (Reversed != PPa->Reversed) ==
                (Pb->Reversed != PPb->Reversed) ? tb->Rank -
                Pb->First->Rank : Pb->Last->Rank - tb->Rank;
            if (nRight < 0)
                nRight += Pb->Size;
        }
        n = nLeft + nMid + nRight + 1;
    }
    return n + 1;
}

#endif
//################################### SegmentSize.c end ###################################

//################################### Sequence.c begin ###################################
#include "Sequence.h"
#include "Segment.h"

/*
 * This file contains the functions FindPermutation and FeasibleKOptMove.
 */

/*  
 * The FindPermutation function finds the permutation p[1:2k] corresponding 
 * to the sequence in which the nodes t[1:2k] occur on the tour.
 *   
 * The nodes are sorted using qsort. The BETWEEN function is used 
 * as comparator.
 *   
 * Postcondition:
 *   
 *     BETWEEN(t[p[i-1]], t[p[i]], t[p[i+1]]) for i = 2, ..., 2k-1
 */

static Node *tp1;

static int Sequence_compare(const void *pa, const void *pb)
{
    return BETWEEN(tp1, t[*(int *) pa], t[*(int *) pb]) ? -1 : 1;
}

void FindPermutation(int k)
{
    int i, j;

    for (i = j = 1; j <= k; i += 2, j++)
        p[j] = SUC(t[i]) == t[i + 1] ? i : i + 1;
    tp1 = t[p[1]];
    qsort(p + 2, k - 1, sizeof(int), Sequence_compare);
    for (j = 2 * k; j >= 2; j -= 2) {
        p[j - 1] = i = p[j / 2];
        p[j] = i & 1 ? i + 1 : i - 1;
    }
    for (i = 1; i <= 2 * k; i++)
        q[p[i]] = i;
}

/*  
 * The FeasibleKOptMove function tests whether the move given by
 * t[1..2k] and incl[1..2k] represents a feasible k-opt move,
 * i.e., making the move on the current tour will result in a tour.
 *   
 * In that case, 1 is returned. Otherwise, 0 is returned. 
 */

int FeasibleKOptMove(int k)
{
    int Count, i;

    FindPermutation(k);
    for (Count = 1, i = 2 * k; (i = q[incl[p[i]]] ^ 1); Count++);
    return Count == k;
}

/*
 * The Cycles function returns the number of cycles that would appear if 
 * the move given by t[1..2k] and incl[1..2k] was made. 
 * In addition, cycle[i] is assigned the number of the cycle that node t[i] 
 * is a part of (an integer from 1 to Cycles).
 */

int Cycles(int k)
{
    int i, j, Count = 0;

    for (i = 1; i <= 2 * k; i++)
        cycle[i] = 0;
    for (i = 1; i <= 2 * k; i++) {
        if (!cycle[p[i]]) {
            Count++;
            j = i;
            do {
                cycle[p[j]] = Count;
                j = q[incl[p[j]]];
                cycle[p[j]] = Count;
                if ((j ^= 1) > 2 * k)
                    j = 1;
            }
            while (j != i);
        }
    }
    return Count;
}

/*
 * The Added function is used to test if an edge, (ta,tb),
 * has been added in the submove under construction.
 */

int Added(const Node * ta, const Node * tb)
{
    return ta->Added1 == tb || ta->Added2 == tb;
}

/* 
 * The Deleted function is used to test if an edge, (ta,tb), 
 * of the tour has been deleted in the submove under construction.
 */

int Deleted(const Node * ta, const Node * tb)
{
    return ta->Deleted1 == tb || ta->Deleted2 == tb;
}

/*
 * The MarkAdded function is used mark an edge, (ta,tb), as added
 * in the submove under construction.
 */

void MarkAdded(Node * ta, Node * tb)
{
    if (!ta->Added1)
        ta->Added1 = tb;
    else if (!ta->Added2)
        ta->Added2 = tb;
    if (!tb->Added1)
        tb->Added1 = ta;
    else if (!tb->Added2)
        tb->Added2 = ta;
}

/*
 * The MarkDeletedfunction is used to mark an edge, (ta,tb), as deleted
 * in the submove under construction.
 */

void MarkDeleted(Node * ta, Node * tb)
{
    if (!ta->Deleted1)
        ta->Deleted1 = tb;
    else if (!ta->Deleted2)
        ta->Deleted2 = tb;
    if (!tb->Deleted1)
        tb->Deleted1 = ta;
    else if (!tb->Deleted2)
        tb->Deleted2 = ta;
}

/*
 * The UnmarkAdded function is used mark the edge, (ta,tb), as not
 * added.
 */

void UnmarkAdded(Node * ta, Node * tb)
{
    if (ta->Added1 == tb)
        ta->Added1 = 0;
    else if (ta->Added2 == tb)
        ta->Added2 = 0;
    if (tb->Added1 == ta)
        tb->Added1 = 0;
    else if (tb->Added2 == ta)
        tb->Added2 = 0;
}

/*
 * The UnmarkDeleted function is used mark the edge, (ta,tb), as not
 * deleted.
 */

void UnmarkDeleted(Node * ta, Node * tb)
{
    if (ta->Deleted1 == tb)
        ta->Deleted1 = 0;
    else if (ta->Deleted2 == tb)
        ta->Deleted2 = 0;
    if (tb->Deleted1 == ta)
        tb->Deleted1 = 0;
    else if (tb->Deleted2 == ta)
        tb->Deleted2 = 0;
}
//################################### Sequence.c end ###################################

//################################### Statistics.c begin ###################################
#include "LKH.h"

static int TrialsMin, TrialsMax, TrialSum, Successes, Updates;
static GainType CostMin, CostMax, CostSum;
static double TimeMin, TimeMax, TimeSum;

void InitializeStatistics()
{
    TrialSum = Successes = Updates = 0;
    CostSum = 0;
    TimeSum = 0.0;
    TrialsMin = INT_MAX;
    TrialsMax = 0;
    TimeMin = DBL_MAX;
    TimeMax = 0;
    CostMin = PLUS_INFINITY;
    CostMax = MINUS_INFINITY;
}

void UpdateStatistics(GainType Cost, double Time)
{
    if (Trial < TrialsMin)
        TrialsMin = Trial;
    if (Trial > TrialsMax)
        TrialsMax = Trial;
    TrialSum += Trial;
    if (Cost <= Optimum)
        Successes++;
    if (Cost < CostMin)
        CostMin = Cost;
    if (Cost > CostMax)
        CostMax = Cost;
    CostSum += Cost;
    if (Time < TimeMin)
        TimeMin = Time;
    if (Time > TimeMax)
        TimeMax = Time;
    TimeSum += Time;
    Updates++;
}

void PrintStatistics()
{
    int _Runs = Updates, _TrialsMin = TrialsMin;
    double _TimeMin = TimeMin;
    GainType _Optimum = Optimum;

    printff("Successes/Runs = %d/%d\n", Successes, Runs);
    if (_Runs == 0)
        _Runs = 1;
    if (_TrialsMin > TrialsMax)
        _TrialsMin = 0;
    if (_TimeMin > TimeMax)
        _TimeMin = 0;
    if (CostMin <= CostMax && CostMin != PLUS_INFINITY) {
        printff
            ("Cost.min = " GainFormat ", Cost.avg = %0.2f, Cost.max = "
             GainFormat "\n", CostMin, (double) CostSum / _Runs, CostMax);
        if (_Optimum == MINUS_INFINITY)
            _Optimum = BestCost;
        if (_Optimum != 0)
            printff
                ("Gap.min = %0.4f%%, Gap.avg = %0.4f%%, Gap.max = %0.4f%%\n",
                 100.0 * (CostMin - _Optimum) / _Optimum,
                 100.0 * ((double) CostSum / _Runs - _Optimum) / _Optimum,
                 100.0 * (CostMax - _Optimum) / _Optimum);

    }
    printff("Trials.min = %d, Trials.avg = %0.1f, Trials.max = %d\n",
            _TrialsMin, 1.0 * TrialSum / _Runs, TrialsMax);
    printff
        ("Time.min = %0.2f sec., Time.avg = %0.2f sec., "
         "Time.max = %0.2f sec.\n",
         fabs(_TimeMin), fabs(TimeSum) / _Runs, fabs(TimeMax));
    printff("Time.total = %0.2f sec.\n", GetTime() - StartTime);
}
//################################### Statistics.c end ###################################

//################################### StoreTour.c begin ###################################
#include "Segment.h"
#include "LKH.h"

/* 
 * The StoreTour function is called each time the tour has been improved by 
 * the LinKernighan function.
 *
 * The function "activates" all nodes involved in the current sequence of moves.
 *
 * It sets OldPred to Pred and OldSuc to Suc for each of these nodes. In this
 * way it can always be determined whether an edge belongs to current starting
 * tour. This is used by the BestMove function to determine whether an edge is
 * excludable.
 *
 * Finally, for each of these nodes the function updates their Cost field.
 * The Cost field contains for each node its minimum cost of candidate edges 
 * not on the tour. The value is used by the BestMove function to decide 
 * whether a tentative non-gainful move should be considered. 
 */

void StoreTour()
{
    Node *t, *u;
    Candidate *Nt;
    int i;

    while (Swaps > 0) {
        Swaps--;
        for (i = 1; i <= 4; i++) {
            t = i == 1 ? SwapStack[Swaps].t1 :
                i == 2 ? SwapStack[Swaps].t2 :
                i == 3 ? SwapStack[Swaps].t3 : SwapStack[Swaps].t4;
            Activate(t);
            t->OldPred = t->Pred;
            t->OldSuc = t->Suc;
            t->OldPredExcluded = t->OldSucExcluded = 0;
            t->Cost = INT_MAX;
            for (Nt = t->CandidateSet; (u = Nt->To); Nt++)
                if (u != t->Pred && u != t->Suc && Nt->Cost < t->Cost)
                    t->Cost = Nt->Cost;
        }
    }
}
//################################### StoreTour.c end ###################################

//################################### SymmetrizeCandidateSet.c begin ###################################
#include "LKH.h"

/* 
 * The SymmetrizeCandidateSet function complements the candidate set such 
 * that every candidate edge is associated with both its two end nodes. 
*/

void SymmetrizeCandidateSet()
{
    Node *From, *To;
    Candidate *NFrom;

    From = FirstNode;
    do {
        for (NFrom = From->CandidateSet; NFrom && (To = NFrom->To);
             NFrom++)
            AddCandidate(To, From, NFrom->Cost, NFrom->Alpha);
    }
    while ((From = From->Suc) != FirstNode);
    ResetCandidateSet();
}
//################################### SymmetrizeCandidateSet.c end ###################################

//################################### TrimCandidateSet.c begin ###################################
#include "LKH.h"

/*
 * The TrimCandidateSet function takes care that each node has 
 * associated at most MaxCandidates candidate edges.                         
 */

void TrimCandidateSet(int MaxCandidates)
{
    Node *From;
    Candidate *NFrom;
    int Count;

    From = FirstNode;
    do {
        Count = 0;
        for (NFrom = From->CandidateSet; NFrom && NFrom->To; NFrom++)
            Count++;
        if (Count > MaxCandidates) {
            From->CandidateSet =
               (Candidate *) realloc(From->CandidateSet,
                                     (MaxCandidates + 1) * sizeof(Candidate));
            From->CandidateSet[MaxCandidates].To = 0;
        }
    } while ((From = From->Suc) != FirstNode);
}
//################################### TrimCandidateSet.c end ###################################

