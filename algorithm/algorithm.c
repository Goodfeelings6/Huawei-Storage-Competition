#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "algorithm.h"
#include "LKHInterface.h"
#include "MaxMatchingByHarryZHR.h" 
#define USING_REAL_TIME_COST  // 是否启用实时的方式算代价，未定义则使用邻接矩阵方式

#define DEBUG_LEVEL 0
// #define DEBUG_LEVEL 1
// #define DEBUG_LEVEL 2

#define INF INT_MAX/200

double MyGetTime(){ // 返回实际时间：秒
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}
int GetObjectValue(const HeadInfo *start, const HeadInfo *end){ /* 返回目标函数值(距离) */
    // return SeekTimeCalculate(start, end);
    return SeekTimeCalculate(start, end)+BeltWearTimes(start, end, NULL);
    // return SeekTimeCalculate(start, end)+BeltWearTimes(start, end, NULL)+MotorWearTimes(start, end);
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
        p->MaxCandidates = 15;
    }
    else if(matDimension <= 5002){ // 5000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        p->POPMUSIC_SampleSize = 20;
        p->POPMUSIC_Solutions = 25;
        p->MaxCandidates = 15;
    }
    else{ // 10000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        p->POPMUSIC_SampleSize = 30;
        p->POPMUSIC_Solutions = 20;
        p->MaxCandidates = 6;
    }
    p->Runs = 1;
    p->TraceLevel = 1;
    p->TimeLimit = DBL_MAX; // 由总时间、一定时间跨度内的改进值共同控制退出即可
    p->TotalTimeLimit = 100; // 最大允许运行时间
    p->ScheduleScoreInSecond = 20;
    p->MoveType = 3;
    p->TimeSpan = 2;
}

// LKH初始解实时算代价
int getCost(const InputParam *input, int len, int i, int j){
    if(i == j)
        return 0;
    else if(i==len-1){ // 磁头节点
        if(j==len-2){ // 磁头节点到虚拟节点的代价为极大值
            return INF;
        }
        else{ // 磁头节点到其他节点的代价为对应寻址时间
            HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            return GetObjectValue(&start, &end);
        }
    }
    else if(i==len-2){ // 虚拟节点
        if(j==len-1){ // 虚拟节点到磁头节点的代价为0
            return 0;
        }
        else{ // 虚拟节点到其他节点的代价为极大值
            return INF;
        }
    }
    else{ // 其他节点
        if(j==len-1){// 其他节点到磁头节点的代价为极大值
            return INF;
        }
        else if(j==len-2){// 其他节点到虚拟节点的代价为0
            return 0;
        }
        else{
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            return GetObjectValue(&start, &end);
        }
    }
}
// 拼接未优化 MaxMatchingByHarryZHR 求解函数
int MaxMatchingByHarryZHRsolver(const InputParam *input, OutputParam *output) {
    /* 生成邻接矩阵 */
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    if (!adjMat) {
        fprintf(stderr, "Memory allocation for adjMat failed\n");
        return -1;
    }
    for (uint32_t i = 0; i < len; ++i) {
        adjMat[i] = (int *)malloc(sizeof(int) * len);
        if (!adjMat[i]) {
            fprintf(stderr, "Memory allocation for adjMat[%d] failed\n", i);
            return -1;
        }
    }

    for (uint32_t i = 0; i < len - 2; ++i) {
        for (uint32_t j = 0; j < len - 2; ++j) {
            if (i == j)
                adjMat[i][j] = INF;
            else {
                HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                adjMat[i][j] = SeekTimeCalculate(&start, &end);
            }
        }
    }

    /* len-2 行 len-2 列是虚拟节点，len-1 行 len-1 列是磁头节点 */
    for (uint32_t i = 0; i < len - 2; ++i) {
        adjMat[len-2][i] = INF;
        adjMat[i][len-2] = 0;

        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = SeekTimeCalculate(&start, &end);
        adjMat[i][len-1] = INF;
    }
    adjMat[len-2][len-1] = 0;
    adjMat[len-2][len-2] = INF;
    adjMat[len-1][len-1] = INF;
    adjMat[len-1][len-2] = INF;

    /*for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            printf("%d\t ", adjMat[i][j]);
        }
        printf("\n");
    }*/

    // int maxCost = 0;
    // for (uint32_t i = 0; i < len; ++i) {
    //     for (uint32_t j = 0; j < len; ++j) {
    //         if (adjMat[i][j] > maxCost && adjMat[i][j] != INF) {
    //             maxCost = adjMat[i][j];
    //         }
    //     }
    // }

    int maxCost = INF;
    //printf("maxcost=%d\n", maxCost);

    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            adjMat[i][j] = maxCost - adjMat[i][j];
           // printf("%d\t", adjMat[i][j]);
        }
       // printf("\n");
    }

// printf("\n");printf("\n");
//     for (uint32_t i = 0; i < len; ++i) {
//         for (uint32_t j = 0; j < len; ++j) {
            
//             printf("%d\t", getCost(input,len,i,j));
//         }
//         printf("\n");
//     }

    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);
    Arr* resultt = solveMaxMatchingByHarryZHR(&matcher, adjMat);
    int *result = (int*)calloc(len, sizeof(int)); ;
    for (uint32_t i = 0; i < len; ++i) 
    { 
        result[i]=getArr(resultt, i);
        //printf("result[i]: %d ", result[i]);
    }    
    int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点
    for(int i=0;i<len;++i)visited[i]=0;
    if (!visited) {
        fprintf(stderr, "Memory allocation for visited failed\n");
        return -1;
    }
    //printf("visited[0]: %d\n", visited[0]);
    //printf("形成的环路：\n");
    int* final_result= (int*)calloc(len, sizeof(int));
    for(int i=0;i<len;++i)final_result[i]=0;
    int count_result=0;
     int loops=0;
    for (int i = len-1; i >= 0; --i) {
        if (!visited[i]) {
            loops++;
            int current = i;
            //printf("环路: %d", current);
            if(current!=len-1&&current!=len-2)
                {   
                    final_result[count_result]=current+1;
                    count_result++;
                }
            while (result[current] != -1 && !visited[result[current]]) {
                visited[current] = 1;
                current = result[current];
                if(current!=len-1&&current!=len-2)
                {   
                    final_result[count_result]=current+1;
                    count_result++;
                }
                //printf(" -> %d", current);
            }
         
            visited[current] = 1;
            //printf("\n");
          
        }
    }
    //for(int i=0;i<10;i++)printf(" %d ", final_result[i]);
    free(visited);
 printf("环路个数: %d\n", loops);
    for (uint32_t i = 0; i < len - 2; ++i) {
        output->sequence[i] = final_result[i] ;
    }
    /*printf("最大匹配算法得到的 io 序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");*/

    // 释放 matcher 和 adjMat 的内存
    freeMaxMatchingByHarryZHR(&matcher);
    for (uint32_t i = 0; i < len; ++i) {
        free(adjMat[i]);
    }
    free(adjMat);

    return 0;
}

//拼接贪心优化 MaxMatchingByHarryZHR 求解函数
int MaxMatchingByHarryZHRsolver2(const InputParam *input, OutputParam *output) {
    /* 生成邻接矩阵 */
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    if (!adjMat) {
        fprintf(stderr, "Memory allocation for adjMat failed\n");
        return -1;
    }
    for (uint32_t i = 0; i < len; ++i) {
        adjMat[i] = (int *)malloc(sizeof(int) * len);
        if (!adjMat[i]) {
            fprintf(stderr, "Memory allocation for adjMat[%d] failed\n", i);
            return -1;
        }
    }

    for (uint32_t i = 0; i < len - 2; ++i) {
        for (uint32_t j = 0; j < len - 2; ++j) {
            if (i == j)
                adjMat[i][j] = INF;
            else {
                HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                adjMat[i][j] = SeekTimeCalculate(&start, &end);
            }
        }
    }

    /* len-2 行 len-2 列是虚拟节点，len-1 行 len-1 列是磁头节点 */
    for (uint32_t i = 0; i < len - 2; ++i) {
        adjMat[len-2][i] = INF;
        adjMat[i][len-2] = 0;

        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = SeekTimeCalculate(&start, &end);
        adjMat[i][len-1] = INF;
    }
    adjMat[len-2][len-1] = 0;
    adjMat[len-2][len-2] = INF;
    adjMat[len-1][len-1] = INF;
    adjMat[len-1][len-2] = INF;

    /*for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) {
            printf("%d\t ", adjMat[i][j]);
        }
        printf("\n");
    }*/

    int maxCost = 0;
    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            if (adjMat[i][j] > maxCost && adjMat[i][j] != INF) {
                maxCost = adjMat[i][j];
            }
        }
    }
    //printf("maxcost=%d\n", maxCost);

    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            adjMat[i][j] = maxCost - adjMat[i][j];
        }
    }

    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);
    Arr* resultt = solveMaxMatchingByHarryZHR(&matcher, adjMat);
    int *result = (int*)calloc(len, sizeof(int)); ;
    for (uint32_t i = 0; i < len; ++i) 
    { 
        result[i]=getArr(resultt, i);
        //printf("result[i]: %d ", result[i]);
    }    
    int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点
    for(int i=0;i<len;++i)visited[i]=0;
    if (!visited) {
        fprintf(stderr, "Memory allocation for visited failed\n");
        return -1;
    }
    

    //printf("visited[0]: %d\n", visited[0]);
    //printf("形成的环路：\n");
    int* final_result= (int*)calloc(len, sizeof(int));
    for(int i=0;i<len;++i)final_result[i]=0;
    int count_result=0;
    int loops=0;
    int current=len-1;//i=len-1,确保从磁头开始构建序列
   while( current!=-1 && visited[current] != 1) {//i=len-1,确保从磁头开始构建序列
       // if (!visited[i]) {
            //int current = i;
            visited[current] = 1;//标记访问
            loops++;
            if(current!=len-1&&current!=len-2)//若当前节点不是磁头节点或者虚拟节点，则放入结果序列中
                {   
                    final_result[count_result]=current+1;
                    count_result++;
                }
            while (result[current] != -1 && !visited[result[current]]) {//继续基于当前节点向后搜索，直到把同一回路（即本回路）的节点全部放入结果序列为止
                visited[current] = 1;//标记访问
                current = result[current];//传递下标索引
                if(current!=len-1&&current!=len-2)//放入结果序列
                {   
                    final_result[count_result]=current+1;
                    count_result++;
                }
              
           }
            
            visited[current] = 1;
            //放完本回路节点之后，遍历所有剩余节点，找到连接路径最短的其他回路节点(这里最短要考虑删边和加边)
            int mincost=INF;
            int nextnode=-1;//记录连接到的下一个回路的节点
            int outtimes=0;
            int intimes=0;
            for (int i = len-1; i >= 0; --i) {
                outtimes++;
                if (!visited[i] ) {
                    intimes++;
                    int j=result[i];//i的下一个节点j才是要判断是否连接的点
                     //printf("mincost=%d\n",mincost);
                    if(mincost>adjMat[i][j]-adjMat[current][j]){
                        mincost=adjMat[i][j]-adjMat[current][j];
                        nextnode=j;
                    }
                }
            }
            //printf("outtimes: %d\n", outtimes);
            //printf("intimes: %d\n", intimes);
            current=nextnode;//把找到的连接路径最短的节点作为下一轮的起始节点
       // }
    }
    //for(int i=0;i<10;i++)printf(" %d ", final_result[i]);
     printf("环路个数: %d\n", loops);


    free(visited);

    for (uint32_t i = 0; i < len - 2; ++i) {
        output->sequence[i] = final_result[i] ;
    }
    // printf("最大匹配算法得到的 io 序列:\n");
    // for (uint32_t i = 0; i < input->ioVec.len; ++i) {
    //     printf("%d ", output->sequence[i]);
    // }
    // printf("\n");

    // 释放 matcher 和 adjMat 的内存
    freeMaxMatchingByHarryZHR(&matcher);
    for (uint32_t i = 0; i < len; ++i) {
        free(adjMat[i]);
    }
    free(adjMat);

    return 0;
}
// LKH 算法
int32_t LKH_Maxmatch(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点


    int *intialTour = (int *)malloc(len * sizeof(int));
    intialTour[0] = len - 1, intialTour[1] = len; 
    int rett=MaxMatchingByHarryZHRsolver(input, output);
    for(int i=2;i<len;i++)intialTour[i]=output->sequence[i-2];
    //for(int i=0;i<len;i++)printf("%d ", intialTour[i]);
    /* 设置固定边, 结点从1开始编号 */
    int fixLen = 1;
    int *fixEdge = (int *)malloc(2 * fixLen * sizeof(int));
    // 虚拟结点到磁头结点的边固定
    fixEdge[0] = len - 1;
    fixEdge[1] = len;

    /* 调用 LKH 求解 */
    /* 确定LKH输入结构体 */
    LKHInput *lkhInput = (LKHInput *)malloc(sizeof(LKHInput));


    lkhInput->adjMat = 0;
    lkhInput->matDimension = len;
    lkhInput->intialTour = intialTour;
    lkhInput->fixEdge = fixEdge;
    lkhInput->fixEdgeLen = fixLen;
    lkhInput->scheduleStartTime = scheduleStartTime;
    lkhInput->lkhParameters = (LKHParameters *)malloc(sizeof(LKHParameters));
    loadDefaultParam(lkhInput->lkhParameters);
    loadUserChangedParam(lkhInput->matDimension, lkhInput->lkhParameters, scheduleStartTime);
#ifdef USING_REAL_TIME_COST
    lkhInput->lkhParameters->OriginInput = input;
#endif
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


    free(intialTour);
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


// LKH 算法
int32_t LKH_Nearest(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点
#ifndef USING_REAL_TIME_COST
    /* 生成邻接矩阵 */
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
                adjMat[i][j] = GetObjectValue(&start, &end);
            }
        }
    }

     /*len-2行len-2列是虚拟节点，len-1行len-1列是磁头节点，为虚拟节点和磁头节点设置代价 */
    for (uint32_t i = 0; i < len - 2; ++i) {
        adjMat[len-2][i] = INF;  // 虚拟节点到其他节点的代价为极大值，则不可能达到其他点
        adjMat[i][len-2] = 0;  // 其他节点到虚拟节点的代价为0
        
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = GetObjectValue(&start, &end);  // 磁头节点到其他节点的代价为对应寻址时间
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
#ifndef USING_REAL_TIME_COST
            int cost = adjMat[current - 1][j - 1];
#endif
#ifdef USING_REAL_TIME_COST
            int cost = getCost(input, len, current - 1, j - 1);
#endif
            if (!visited[j] && cost < min_dist){
                min_dist = cost;
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
#ifndef USING_REAL_TIME_COST
    lkhInput->adjMat = adjMat;
#endif
#ifdef USING_REAL_TIME_COST
    lkhInput->adjMat = 0;
#endif
    lkhInput->matDimension = len;
    lkhInput->intialTour = intialTour;
    lkhInput->fixEdge = fixEdge;
    lkhInput->fixEdgeLen = fixLen;
    lkhInput->scheduleStartTime = scheduleStartTime;
    lkhInput->lkhParameters = (LKHParameters *)malloc(sizeof(LKHParameters));
    loadDefaultParam(lkhInput->lkhParameters);
    loadUserChangedParam(lkhInput->matDimension, lkhInput->lkhParameters, scheduleStartTime);
#ifdef USING_REAL_TIME_COST
    lkhInput->lkhParameters->OriginInput = input;
#endif
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
#ifndef USING_REAL_TIME_COST
    for (int i = 0; i < len; i++) {
        free(adjMat[i]);  // 逐行释放
    }
    free(adjMat);
#endif
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
    uint32_t len = output->len + 1; //增加了一个磁头起始位置节点
    /* 生成邻接矩阵 */
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
                adjMat[i][j] = GetObjectValue(&start, &end);;
            }
        }
    }

     /*len-1行len-1列是磁头节点，为磁头节点设置代价 */
    for (uint32_t i = 0; i < len - 1; ++i) {  
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = GetObjectValue(&start, &end);;  // 磁头节点到其他节点的代价为对应寻址时间
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
            int cost = adjMat[current - 1][j - 1];
            if(!visited[j] && cost < min_dist){
                min_dist = cost;
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
    return LKH_Maxmatch(input, output);
    //return LKH_Nearest(input, output);
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
