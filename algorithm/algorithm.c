#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "algorithm.h"
#include "LKHInterface.h"
#include "MaxMatchingByHarryZHR.h" 
#include"priorityQueue.h"
#define USING_REAL_TIME_COST  // 是否启用实时的方式算代价，未定义则使用邻接矩阵方式

#define DEBUG_LEVEL 0
//#define DEBUG_LEVEL 1
// #define DEBUG_LEVEL 2

#define INF INT_MAX/200

// 权重控制系数
double alpha = 0.3; // 读时延权重
double beta = 0.5; // 带体磨损权重
double gama = 0.2; // 电机磨损权重
double base_totaltime = 0;     // 基线读时延
double base_tapeBeltWear = 0;  // 基线带体磨损
double base_tapeMotorWear = 0; // 基线电机磨损
double maxbase = 0;            // 缩放系数
double final_alpha = 0;
double final_beta = 0;
double final_gama = 0;
typedef struct
{
    uint32_t id;         // IO序号
    uint32_t wrap;       // 起始wrap
    uint32_t startLpos;  // 起始lpos
    uint32_t endLpos;    // 结束lpos
    uint32_t visit;      //是否被选
} IOUintnn;

// 调用SCAN基线并分析以获取权重控制系数
void getbaseline(const InputParam *input, OutputParam *output)
{
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();
    Scan(input, output);
    /* 基线读时延 */
    base_totaltime += MyGetTime() - scheduleStartTime; // 加上SCAN算法排序时间
    AccessTime accessTime = {0};
    TotalAccessTime(input, output, &accessTime);
    base_totaltime += accessTime.addressDuration; // 加上读寻址时间
    base_totaltime += accessTime.readDuration;    // 加上读数据时间
    /* 基线带体磨损 */
    TapeBeltSegWearInfo segWearInfo = {0};
    base_tapeBeltWear = TotalTapeBeltWearTimes(input, output, &segWearInfo);
    /* 基线电机磨损 */
    base_tapeMotorWear = TotalMotorWearTimes(input, output);
    // printf("base_totaltime = %lf\n", base_totaltime);
    // printf("base_tapeBeltWear = %lf\n", base_tapeBeltWear);
    // printf("base_tapeMotorWear = %lf\n", base_tapeMotorWear);

    if (base_tapeBeltWear < base_totaltime){
        maxbase = base_totaltime;
    }
    else{
        maxbase = base_tapeBeltWear;
    }
    maxbase = maxbase*10;
    base_totaltime /= maxbase;
    base_tapeBeltWear /= maxbase;
    base_tapeMotorWear /= maxbase;
    // printf("maxbase = %lf\n", maxbase);
    // printf("new base_totaltime = %lf\n", base_totaltime);
    // printf("new base_tapeBeltWear = %lf\n", base_tapeBeltWear);
    // printf("new base_tapeMotorWear = %lf\n", base_tapeMotorWear);

    /* 检测算例场景 
       备份归档场景backup：alpha=0.3,beta=0.5,gama=0.2 
       高性能场景hdd：alpha=0.5,beta=0.3,gama=0.2 
    */
    int Nosequential_count=0;//计数非顺序IO节点
    for (uint32_t i = 1; i < input->ioVec.len; ++i)
    {
        if (input->ioVec.ioArray[i].wrap < input->ioVec.ioArray[i - 1].wrap) // 如果后一个wrap小于前一个 肯定非顺序
        {
            Nosequential_count++;
            // printf("Nosequential_id = %d\n", input->ioVec.ioArray[i].id);
        }
        else if (input->ioVec.ioArray[i].wrap == input->ioVec.ioArray[i - 1].wrap)
        {
            if (input->ioVec.ioArray[i].wrap % 2 == 0) // 正向
            {
                if (input->ioVec.ioArray[i].startLpos < input->ioVec.ioArray[i - 1].startLpos)
                    Nosequential_count++;
                // printf("Nosequential_id = %d\n", input->ioVec.ioArray[i].id);
            }
            else // 反向
            {
                if (input->ioVec.ioArray[i].startLpos > input->ioVec.ioArray[i - 1].startLpos)
                    Nosequential_count++;
                // printf("Nosequential_id = %d\n", input->ioVec.ioArray[i].id);
            }
        }
        // 是否 hdd 场景
        if (Nosequential_count >= input->ioVec.len * 0.1)
        {
            alpha = 0.5; // 读时延权重
            beta = 0.3;  // 带体磨损权重
            gama = 0.2;  // 电机磨损权重
            break;
        }
    }
    // printf("Nosequential_count = %d\n", Nosequential_count);
    // if(alpha==0.3)
    //     printf("backup\n");
    // else
    //     printf("hdd\n");

    // 计算最终权重
    final_alpha = alpha / base_totaltime;
    final_beta = beta / base_tapeBeltWear;
    final_gama = gama / base_tapeMotorWear;
}

double MyGetTime(){ // 返回实际时间：秒
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

/* 返回目标函数值(距离) */
int GetObjectValue(const HeadInfo *start, const HeadInfo *end){
    //double value = SeekTimeCalculate(start, end) * alpha / base_totaltime + BeltWearTimes(start, end, NULL) * beta / base_tapeBeltWear + MotorWearTimes(start, end) * gama / base_tapeMotorWear;
    double value = SeekTimeCalculate(start, end)*final_alpha+MotorWearTimes(start, end)*final_gama;
    // double value = BeltWearTimes(start, end, NULL) * beta / base_tapeBeltWear +MotorWearTimes(start, end)*gama/base_tapeMotorWear;
    // if(value<100)
    //     printf("value:%f\n", value);
    return (int)value;

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

}
/* 返回目标函数值(距离) */
int GetObjectValue_LKH(const HeadInfo *start, const HeadInfo *end){
    double value = SeekTimeCalculate(start, end) * final_alpha + BeltWearTimes(start, end, NULL) * final_beta + MotorWearTimes(start, end) * final_gama;
    return (int)value;
}
// 设置需要调整的 LKH 的参数
void loadUserChangedParam(int matDimension, LKHParameters *p, double scheduleStartTime){
    if(matDimension <= 102){ // 10,50,100
        // 默认值
    }
    else if(matDimension <= 1002){ // 1000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag1000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 10;
        p->TimeSpan = 0.2;
    }
    else if(matDimension <= 2002){ // 2000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag2000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 5;
        p->TimeSpan = 0.2;
    }
    else if(matDimension <= 5002){ // 5000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag5000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 4;
        p->TimeSpan = 0.3;
    }
    else{ // 10000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag10000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 3;
        p->TimeSpan = 0.5;
    }
    p->Runs = 1;
    p->TraceLevel = 0;
    p->TimeLimit = DBL_MAX; // 由总时间、一定时间跨度内的改进值共同控制退出即可
    p->TotalTimeLimit = 40; // 最大允许运行时间
    p->PenaltyScoreInSecond = alpha*1000/base_totaltime + ((2-(matDimension-2)/10000.0)/200.0)*maxbase; /* 20s后每多算1秒实际罚分 (相对Cost) */
    p->ScheduleScoreInSecond = alpha*1000/base_totaltime + (((matDimension-2)/10000.0)/200.0)*maxbase;/* 20s内每少算1秒实际加分 (相对Cost) */
    // printf("p->ScheduleScoreInSecond = %lf\n", p->ScheduleScoreInSecond);
    // printf("p->PenaltyScoreInSecond = %lf\n", p->PenaltyScoreInSecond);
    p->MoveType = 3;
    
}

// 实时算代价
int getCost(const InputParam *input, int len, int i, int j){
    if(i == j)
        return INF;
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

// 邻接矩阵算代价
void getAdjMat(const InputParam *input, int len, int** adjMat){
    //为了节点能和id号对应，还是将它们从0开始对应行列
    for (uint32_t i = 0; i < len - 2; ++i) {
        for (uint32_t j = 0; j < len - 2; ++j) {
            if(i == j)
                adjMat[i][j] = INF;
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
    adjMat[len-2][len-2] = INF; // 虚拟节点到自身的代价为极大值
    adjMat[len-1][len-1] = INF; // 磁头节点到自身的代价为极大值
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
}

// 邻接表算代价
void getAdjList1(const InputParam* input, Graph* graph) {
    uint32_t len = graph->n;

 for (uint32_t i = 0; i < len - 2; ++i) {
        NodeCost *nodeCosts = (NodeCost *)malloc((len - 1) * sizeof(NodeCost));  // len-1是因为不考虑节点到自己的代价
        int count = 0;
        for (uint32_t j = 0; j < len - 2; ++j) {
            if(i == j)
            continue;
            
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            nodeCosts[count].node=j;
            nodeCosts[count].cost =-GetObjectValue(&start, &end);
            count++;
        }

        // 对代价数组按代价升序排序
        qsort(nodeCosts, count, sizeof(NodeCost), compareEdges);

        // 打印当前节点i的最近十个节点
        // printf("Node %d's nearest 10 nodes: ", i);
        for (int k = 0; k < 500 && k < count; ++k) {
            addEdge(graph, i, nodeCosts[k].node, nodeCosts[k].cost);
           // printf("%d (cost: %d) ", nodeCosts[k].node, nodeCosts[k].cost);
        }
        // printf("\n");

        // 释放内存
        free(nodeCosts);

    }
    // 处理虚拟节点和磁头节点的代价
    // 虚拟节点 (len-2) 到其他节点的代价为 INF

    for (uint32_t i = 0; i < len - 2; ++i) {
        //addEdge(graph, len - 2, i, -INF);  // 虚拟节点到其他节点的代价为 INF
        addEdge(graph, i, len - 2, 0);    // 其他节点到虚拟节点的代价为 0

        // 磁头节点 (len-1) 到其他节点的代价为计算后的值
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        Cost headCost = GetObjectValue(&start, &end);
        addEdge(graph, len - 1, i, -headCost);  // 磁头节点到任务的代价
        //addEdge(graph, i, len - 1, -INF);      // 其他节点到磁头节点的代价为 INF
    }

    // 虚拟节点到磁头节点的代价为 0
    addEdge(graph, len - 2, len - 1, 0);  // 虚拟节点到磁头节点的代价为 0
    //addEdge(graph, len - 1, len - 2, -INF);  // 磁头节点到虚拟节点的代价为 INF

    // // 虚拟节点和磁头节点到自身的代价为 INF
    // addEdge(graph, len - 2, len - 2, -INF);  // 虚拟节点到自己
    // addEdge(graph, len - 1, len - 1, -INF);  // 磁头节点到自己
}

void getAdjList2(const InputParam* input, Graph* graph) {
    uint32_t len = graph->n;

    for (uint32_t i = 0; i < len - 1; ++i) {
        NodeCost heap[HEAP_SIZE];  // 大顶堆数组
        int heapSize = 0;
        int iwrap=input->ioVec.ioArray[i].wrap;
        int index_samewrap_close=0;
        int mindis=INF;
        for (uint32_t j = 0; j < len - 1; ++j) {
            int jwrap=input->ioVec.ioArray[j].wrap;
            if (i == j) {
                //addEdge(graph, i, i, -INF+1);
                continue;
            }
          
            if(iwrap==jwrap){//对同wrap的处理
                if(iwrap%2==0){//正向
                    if(input->ioVec.ioArray[j].startLpos < input->ioVec.ioArray[i].endLpos){
                        continue;
                    }
                    else if(mindis >= input->ioVec.ioArray[j].startLpos - input->ioVec.ioArray[i].endLpos){
                        mindis=input->ioVec.ioArray[j].startLpos - input->ioVec.ioArray[i].endLpos;
                        index_samewrap_close=j;//同wrap里i最近的IO块是j
                    }
                    
                }else{
                    if(input->ioVec.ioArray[j].startLpos > input->ioVec.ioArray[i].endLpos){
                        continue;
                    }
                    else if(mindis >= input->ioVec.ioArray[i].endLpos -input->ioVec.ioArray[j].startLpos){
                        mindis=input->ioVec.ioArray[i].endLpos -input->ioVec.ioArray[j].startLpos;
                        index_samewrap_close=j;//同wrap里i最近的IO块是j
                    }
                }
            continue;
            }
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            int cost = GetObjectValue_LKH(&start, &end);

            NodeCost newElement = {j, cost};
            insertMaxHeap(heap, &heapSize, newElement);
        }

        // 将堆中元素添加到图的邻接表中
        for (int k = 0; k < heapSize; ++k) {
            addEdge(graph, i, heap[k].node, -heap[k].cost);
            //printf("Node: %d, Cost: %lf\t", heap[k].node, heap[k].cost);
        }
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[index_samewrap_close].wrap, input->ioVec.ioArray[index_samewrap_close].startLpos, HEAD_RW};
            int cost = GetObjectValue_LKH(&start, &end);
            addEdge(graph, i, index_samewrap_close, -cost);
    }

    // 处理虚拟节点和磁头节点的代价
    for (uint32_t i = 0; i < len - 1; ++i) {
        addEdge(graph, i, len - 1, 0);

        // HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        // HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        // int headCost = GetObjectValue(&start, &end);
        // addEdge(graph, len - 1, i, -headCost);
    }
    addEdge(graph, len - 1, len - 2, 0);
}
void getAdjList_GC(const InputParam* input, Graph* graph) {
    uint32_t len = graph->n;

    for (uint32_t i = 0; i < len - 1; ++i) {
        NodeCost heap[HEAP_SIZE];  // 大顶堆数组
        int heapSize = 0;
        int iwrap=input->ioVec.ioArray[i].wrap;
        int index_samewrap_close=0;
        int mindis=INF;
        for (uint32_t j = 0; j < len - 1; ++j) {
            int jwrap=input->ioVec.ioArray[j].wrap;
            if (i == j) {
                //addEdge(graph, i, i, -INF+1);
                continue;
            }
          
            if(iwrap==jwrap){//对同wrap的处理
                if(iwrap%2==0){//正向
                    if(input->ioVec.ioArray[j].startLpos < input->ioVec.ioArray[i].endLpos){
                        continue;
                    }
                    else if(mindis >= input->ioVec.ioArray[j].startLpos - input->ioVec.ioArray[i].endLpos){
                        mindis=input->ioVec.ioArray[j].startLpos - input->ioVec.ioArray[i].endLpos;
                        index_samewrap_close=j;//同wrap里i最近的IO块是j
                    }
                    
                }else{
                    if(input->ioVec.ioArray[j].startLpos > input->ioVec.ioArray[i].endLpos){
                        continue;
                    }
                    else if(mindis >= input->ioVec.ioArray[i].endLpos -input->ioVec.ioArray[j].startLpos){
                        mindis=input->ioVec.ioArray[i].endLpos -input->ioVec.ioArray[j].startLpos;
                        index_samewrap_close=j;//同wrap里i最近的IO块是j
                    }
                }
            continue;
            }
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            int cost = GetObjectValue(&start, &end);

            NodeCost newElement = {j, cost};
            insertMaxHeap(heap, &heapSize, newElement);
        }

        // 将堆中元素添加到图的邻接表中
        for (int k = 0; k < heapSize; ++k) {
            addEdge(graph, i, heap[k].node, -heap[k].cost);
            //printf("Node: %d, Cost: %lf\t", heap[k].node, heap[k].cost);
        }
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[index_samewrap_close].wrap, input->ioVec.ioArray[index_samewrap_close].startLpos, HEAD_RW};
            int cost = GetObjectValue(&start, &end);
            addEdge(graph, i, index_samewrap_close, -cost);
    }

    // 处理虚拟节点和磁头节点的代价
    for (uint32_t i = 0; i < len - 1; ++i) {
        addEdge(graph, i, len - 1, 0);

        // HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        // HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        // int headCost = GetObjectValue(&start, &end);
        // addEdge(graph, len - 1, i, -headCost);
    }
    addEdge(graph, len - 1, len - 2, 0);
}

// 贪心选边策略
int Greedy_chooseEdge(const InputParam *input, OutputParam *output) {
    
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
    
    // 初始化图
    GraphGC graph;
    initGraphGC(&graph, len);
    NodeCost* queues[len];  // 每个点一个优先队列
    int queueSizes[len];    // 每个队列的当前大小
    int queuehead[len];    // 每个队列当前队首
    // 初始化优先队列
    for (int i = 0; i < len; i++) {
        queues[i] = (NodeCost*)malloc(sizeof(NodeCost) * PriorityQueuesize); // 每个队列大小最多为 k
        queueSizes[i] = 0;                                  // 初始队列为空
        queuehead[i] = 0;
    }

    // 构建每个点构建优先队列
    for(int i=0;i<len-1;i++){
        for(int j=0;j<len-1;j++){
            HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            int cost = GetObjectValue(&start, &end);

            NodeCost newElement = {j, cost};
            insertMinPriorityQueue(queues[i], &queueSizes[i], newElement);
        }
    }
   
    for (uint32_t i = 0; i < len - 2; ++i) {
        NodeCost newElement = {len-1, 0};
        insertMinPriorityQueue(queues[i], &queueSizes[i],  newElement); // 其他节点到虚拟节点的代价为0
        NodeCost newElement2 = {i,INF};
        insertMinPriorityQueue(queues[len-1], &queueSizes[len-1], newElement2);// 虚拟节点到其他节点的代价为INF

        // 磁头节点 (len-2) 到其他节点的代价为计算后的值
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        int headCost = GetObjectValue(&start, &end);
        NodeCost newElement3 = {i, headCost};
        insertMinPriorityQueue(queues[len-2], &queueSizes[len-2], newElement3); // 其他节点到虚拟节点的代价为0
      
    }
    NodeCost newElement = {len-2,0};
    // 虚拟节点到磁头节点的代价为 0
    insertMinPriorityQueue(queues[len-1], &queueSizes[len-1],newElement); // 虚拟节点到磁头节点的代价为0  

    printf("优先队列构建完成\n");
    
    int* outvisited = (int*)calloc(len, sizeof(int)); // 每个节点出度
    int* invisited = (int*)calloc(len, sizeof(int)); // 每个节点入度
    
    int current=1;
    for(int i=0;i<len-1;i++){//做len-1次选边动作，找出len-1条边连接len个结点构成路径
        int minedge=INF;
        int index=-1;
        printf("%d\n",i);
        for(int j=len-1;j>=0;j--){//遍历所有节点的相连边，找没有被选入的最大边
            printf("queuehead[%d]=%d\n",j,queuehead[j]);
        // 遍历每个队列头，避免越界访问
            printf("outvisited[%d]=%d,\n",j,outvisited[j]);
            if(outvisited[j]==0 && invisited[queues[j][queuehead[j]].node]==0){//连出边不能有出度，连线边不能有入度
           
                if( current < len-1 &&(invisited[j]==0 || outvisited[queues[j][queuehead[j]].node]==0)){//连出边不能有入度，或者连线边不能有出度，禁止成环
                    if(minedge>queues[j][queuehead[j]].cost){
                      minedge=queues[j][queuehead[j]].cost;
                      index=j;
                    }
                }else if( current == len-1){
                    if(minedge>queues[j][queuehead[j]].cost){
                      minedge=queues[j][queuehead[j]].cost;
                      index=j;
                    }
                }
            }
        }
        printf("%d\n",index);
        printf("选入的边为%d->%d,cost=%d\n",index, queues[index][queuehead[index]].node, queues[index][queuehead[index]].cost);
        addEdgeGC(&graph, index, queues[index][queuehead[index]].node);//将边选入图中
        outvisited[index]++;//增加连出边出度
        invisited[queues[index][queuehead[index]].node]++;//增加连入边入度
        queuehead[index]++;
    }
    
    printf("输出序列：\n");
    EdgeGC *e = graph.lists[len-2].head;
        while (e) {
            EdgeGC *tmp = e;
            e = e->next;
            output->sequence[current++]=e->id;
            printf("%d\t",e->id);
            free(tmp);
    }
    
    return 0;
}

// LKH 算法 + 最大匹配(MM)生成初始解
int32_t LKH_MM(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点

    int *intialTour = (int *)malloc(len * sizeof(int));
    intialTour[0] = len - 1, intialTour[1] = len;
    MaxMatching_Greedy(input, output);
    // MaxMatching_Greedy(input, output);
    for (int i = 2; i < len; i++)
        intialTour[i] = output->sequence[i - 2];
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

    return ret;
}
// LKH 算法 + 扫描构造生成初始解
int32_t LKH_nn1(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点

    int *intialTour = (int *)malloc(len * sizeof(int));
    intialTour[0] = len - 1, intialTour[1] = len;
    NearestNeighbor1(input, output);
    // MaxMatching_Greedy(input, output);
    for (int i = 2; i < len; i++)
        intialTour[i] = output->sequence[i - 2];
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

    return ret;
}
int32_t LKH_nn2(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点

    int *intialTour = (int *)malloc(len * sizeof(int));
    intialTour[0] = len - 1, intialTour[1] = len;
    NearestNeighbor2(input, output);
    // MaxMatching_Greedy(input, output);
    for (int i = 2; i < len; i++)
        intialTour[i] = output->sequence[i - 2];
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

    return ret;
}

// LKH 算法 + 最近邻(NN)生成初始解
int32_t LKH_NN(const InputParam *input, OutputParam *output)
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
    getAdjMat(input, len, adjMat);
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
            if(!visited[j]){
#ifndef USING_REAL_TIME_COST
                int cost = adjMat[current - 1][j - 1];
#else
                int cost = getCost(input, len, current - 1, j - 1);
#endif
                if (!visited[j] && cost < min_dist){
                    min_dist = cost;
                    node = j;
                }
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
#else
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

    return ret;
}

// LKH 算法 + 贪心插入(GI)生成初始解
int32_t LKH_GI(const InputParam *input, OutputParam *output)
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
    getAdjMat(input, len, adjMat);
#endif

    /* 基于贪心插入构造初始解， 注意巡回路径从1开始编号*/
    int *Next = (int *)calloc(len + 1, sizeof(int)); //记录当前城市的下一个城市
    // 从虚拟结点出发, 再到磁头结点
    Next[len - 1] = len;
    for (int i = 1; i <= len-2; ++i){ // 每次找一个城市
        int minDeltaCost = INF;
        int insertPrev = -1;
        // 评估磁头往后的k个插入位置
        int prev = len; // 从磁头开始
        while(prev != 0){
            int deltaCost;
            if(Next[prev] == 0) // 末尾位置插入
#ifndef USING_REAL_TIME_COST
                deltaCost = adjMat[prev-1][i-1];
#else
                deltaCost = getCost(input, len, prev-1, i-1);
#endif
            else // 中间位置插入
#ifndef USING_REAL_TIME_COST
                deltaCost = adjMat[prev-1][i-1] + adjMat[i-1][Next[prev]-1] - adjMat[prev-1][Next[prev]-1];
#else
                deltaCost = getCost(input, len, prev-1, i-1) + getCost(input, len, i-1, Next[prev]-1) - getCost(input, len, prev-1, Next[prev]-1);
#endif
            if (deltaCost < minDeltaCost){
                minDeltaCost = deltaCost;
                insertPrev = prev;
            }
            prev = Next[prev];
        }
        // 当前城市插入得到的最优位置
        if(Next[insertPrev] == 0){ // 末尾位置插入
            Next[insertPrev] = i;
        }
        else{ // 中间位置插入
            Next[i] = Next[insertPrev];
            Next[insertPrev] = i;
        }
    }
    // 路径链表拷贝至初始解数组
    int *intialTour = (int *)malloc(len * sizeof(int));
    int curr = len-1, idx = 0;
    while(curr != 0){
        intialTour[idx++] = curr;
        // printf("%d ", curr);
        curr = Next[curr];
    }
    // printf("\n");

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
    lkhInput->lkhParameters->OriginInput = input;
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
    free(Next);
    free(intialTour);
    free(fixEdge);
    free(lkhInput->lkhParameters);
    free(lkhInput);
    free(lkhOutput->tourResult);
    free(lkhOutput);

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
    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点
#ifndef USING_REAL_TIME_COST
    /* 生成邻接矩阵 */
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    for (uint32_t i = 0; i < len; ++i){
        adjMat[i] = (int *)malloc(sizeof(int) * len);
    }
    getAdjMat(input, len, adjMat);
#endif

    /* 贪心构造, 注意巡回路径从1开始编号 */
    int *tour = (int *)malloc(len * sizeof(int));
    int* visited = (int*) calloc(len + 1, sizeof(int)); // 标记已经到达过的城市，数组初始化为0
    // 从虚拟结点出发, 再到磁头结点
    tour[0] = len - 1, tour[1] = len; 
    visited[len-1] = 1, visited[len] = 1;
    int current = len; // 磁头
    for (int i = 2; i < len; ++i){
        int min_dist = INF;
        int node = -1;
        // 找到最近的未访问城市
        for (int j = 1; j <= len; ++j){
            if(!visited[j]){
#ifndef USING_REAL_TIME_COST
                int cost = adjMat[current - 1][j - 1];
#else
                int cost = getCost(input, len, current - 1, j - 1);
#endif
                if (cost < min_dist){
                    min_dist = cost;
                    node = j;
                }
            }
        }
        // 更新当前城市
        visited[node] = 1;
        tour[i] = node;
        current = node;
    }
    
    // 赋值到输出数据结构
    for (int i = 0; i < len - 2; ++i){ 
        output->sequence[i] = tour[i+2];
    }
#if DEBUG_LEVEL >= 1
    printf("最近邻算法得到的io序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
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
    free(tour);
    free(visited);

    return 0;
}
int findFirstGreaterOrEqual(IOUintnn * arr, int n, int target) {
    int left = 0, right = n - 1;
    int result = -1; // 初始化结果为 -1，表示未找到
    //printf("result: %d\n",result);
    //printf("target: %d\n",target);
    while (left <= right) {
        int mid = left + (right - left) / 2;
        //printf("target: %d  arr[mid].startLpos:%d  \n",target,arr[mid].startLpos);
        if (arr[mid].startLpos > target) {
            //printf("第一个大于的数 %d\n",arr[mid].startLpos);
            result = mid; // 更新结果为当前 mid
            right = mid - 1; // 继续向左找，确保是第一个
        } else {
            left = mid + 1; // 向右找更大的值
        }
    }
    //printf("result: %d\n",result);
    return result; // 返回第一个满足条件的下标（或 -1 表示未找到）
}

int findFirstLessOrEqual(IOUintnn * arr, int n, int target) {
    int left = 0, right = n - 1;
    int result = -1; // 初始化结果为 -1，表示未找到

    while (left <= right) {
        int mid = left + (right - left) / 2;

        if (arr[mid].startLpos < target) {
            //printf("第%d个小于%d的数 %d\n",mid,target,arr[mid].startLpos);
            result = mid; // 更新结果为当前 mid
            left = mid + 1; // 继续向右找，确保是第一个
        } else {
            right = mid - 1; // 向左找更小的值
        }
    }

    return result; // 返回第一个满足条件的下标（或 -1 表示未找到）
}

int nn_compare(const void *a, const void *b) {
    IOUintnn *ioA = (IOUintnn *)a;
    IOUintnn *ioB = (IOUintnn *)b;

    if (ioA->startLpos < ioB->startLpos) {
        return -1; // ioA 在 ioB 之前
    } else if (ioA->startLpos >= ioB->startLpos) {
        return 1;  // ioA 在 ioB 之后
    } else {
        return 0;  // 两者相等
    }
}

void deleteElement(IOUintnn * arr, int *n, int index) {
    if (index < 0 || index >= *n) {
        //printf("Invalid index to delete.\n");
        return;
    }

    // 将后续元素左移覆盖
    for (int i = index; i < *n - 1; i++) {
        arr[i] = arr[i + 1];
    }

    // 减少数组长度
    (*n)--;
}

int nonNegative(int num) {
    return num >= 0 ? num : 0;
}
int32_t NearestNeighbor1(const InputParam *input, OutputParam *output){    
    // 创建 ioArray 的副本，用于排序
    //IOUintnn *nn_IOArray = (IOUintnn *)malloc(input->ioVec.len * sizeof(IOUintnn));
    int len1=0;
    int len0=0;
    for (uint32_t i = 0; i < input->ioVec.len; ++i) 
    {
        if(input->ioVec.ioArray[i].wrap%2==1)len1++;
        else len0++;
    }
    IOUintnn *nnodd_IOArray = (IOUintnn *)malloc(len1* sizeof(IOUintnn));;
    IOUintnn *nneven_IOArray  = (IOUintnn *)malloc(len0 * sizeof(IOUintnn));;

    uint32_t oddIndex = 0;
    uint32_t evenIndex = 0;
    
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        if(input->ioVec.ioArray[i].wrap%2==1)
        {
            nnodd_IOArray[oddIndex].endLpos = input->ioVec.ioArray[i].endLpos+5;
            nnodd_IOArray[oddIndex].startLpos = input->ioVec.ioArray[i].startLpos+5;
            nnodd_IOArray[oddIndex].wrap=input->ioVec.ioArray[i].wrap;
            nnodd_IOArray[oddIndex].id=input->ioVec.ioArray[i].id;
            nnodd_IOArray[oddIndex].visit=0;
            oddIndex++;
        }
        else
        {
            nneven_IOArray[evenIndex].endLpos = input->ioVec.ioArray[i].endLpos+5;
            nneven_IOArray[evenIndex].startLpos = input->ioVec.ioArray[i].startLpos+5;
            nneven_IOArray[evenIndex].wrap=input->ioVec.ioArray[i].wrap;
            nneven_IOArray[evenIndex].id=input->ioVec.ioArray[i].id;
            nneven_IOArray[evenIndex].visit=0;
            evenIndex++;
        }
    }
    
    // 使用 qsort 对副本进行排序
    qsort(nnodd_IOArray, oddIndex, sizeof(IOUintnn), nn_compare);
    qsort(nneven_IOArray, evenIndex, sizeof(IOUintnn), nn_compare);
    //input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status;
    //for (uint32_t i = 0; i < len0; ++i)
    //printf("偶数第%d个节点%dstarlops: %d\n",i,nneven_IOArray[i].id,nneven_IOArray[i].startLpos);
    //for (uint32_t i = 0; i < len1; ++i)
    //printf("奇数第%d个节点%dstarlops: %d\n",i,nnodd_IOArray[i].id,nnodd_IOArray[i].startLpos);
    uint32_t seqIndex = 0;
    double target=input->headInfo.lpos;
    int xx=1;
    //printf("循环2:\n");
    if(input->headInfo.wrap%2==0)
    {
        xx=0;
    }
    /*
    if(input->headInfo.wrap%2==0)
    {   printf("算法进入了0:\n");
        int index = findFirstGreaterOrEqual(nneven_IOArray, evenIndex, input->headInfo.lpos);
        printf("index: %d\n",index);
        nneven_IOArray[index].visit=1;
        target=nneven_IOArray[index].endLpos;
        output->sequence[seqIndex++]=nneven_IOArray[index].id;
        xx=0;
    }
    else
    {   printf("算法进入了1:\n");
        int index = findFirstLessOrEqual(nnodd_IOArray, oddIndex, input->headInfo.lpos);
        nnodd_IOArray[index].visit=1;
        target=nnodd_IOArray[index].endLpos;
        output->sequence[seqIndex++]=nnodd_IOArray[index].id;
        xx=1;
    }
    */
    //printf("第一个: %d\n",output->sequence[0]);
    //printf("seqIndex: %d\n",seqIndex);
    //printf("target: %d\n",target);
    /**/
    int evenIndextag=1;
    int oddIndextag=1;
    int index;
    while(seqIndex<input->ioVec.len)
    //for(int i=0;i<10;i++)
    {   //printf("循环 %d xx:%d:\n",seqIndex,xx);
        //printf("xx %d:\n",xx);
        //printf("evenIndex %d:\n",evenIndex);
        //printf("evenIndex: %d || xx :%d\n",evenIndex,xx);
        if(xx==0)
        {   //printf("evenIndex: %d || xx :%d\n",evenIndex,xx);
            if(evenIndextag!=0)
            {
            index = findFirstGreaterOrEqual(nneven_IOArray, evenIndex, target);
            //printf("index: %d || xx :%d\n",index,xx);
            if(index!=-1)
            {   //printf("算法进入了1:\n");
                nneven_IOArray[index].visit=1;
                target=nneven_IOArray[index].endLpos;
                output->sequence[seqIndex++]=nneven_IOArray[index].id;
                deleteElement(nneven_IOArray, &evenIndex, index);
                len0--;
                
                if(index==evenIndex-1)
                {   //printf("算法进入了3:\n");
                    if(oddIndex==0)target=(nneven_IOArray[0].startLpos)-1;
                    else target=nnodd_IOArray[oddIndex-1].startLpos+1;
                }
               
            }
            else if(index==-1)
            {
                
                //printf("算法进入了2:\n");
                 if(evenIndex==0)evenIndextag=0;
                 //printf("oddIndex: %d\n",oddIndex);
                 if(oddIndex!=0){xx=1;target=nnodd_IOArray[oddIndex-1].startLpos+1;}
                 else 
                 {
                    target=nneven_IOArray[0].startLpos*0.5;;
                 }
            }
             }

        }
        else
        {   
              if(oddIndextag!=0)
             {
            //printf("算法错误:\n");
            index = findFirstLessOrEqual(nnodd_IOArray, oddIndex, target);
             //printf("index: %d || xx :%d\n",index,xx);
            if(index!=-1)
            {   //printf("算法进入了1:\n");
                nnodd_IOArray[index].visit=1;
                target=nnodd_IOArray[index].endLpos;
                output->sequence[seqIndex++]=nnodd_IOArray[index].id;
                deleteElement(nnodd_IOArray, &oddIndex, index);
                len1--;
                if(index==0)
                {   //printf("算法进入了2%d:\n,oddIndex");
                    if(evenIndex==0){target=nnodd_IOArray[nonNegative(oddIndex-1)].startLpos+1;}
                    else {target=nneven_IOArray[index].startLpos-1;;}
                    //printf("target: %d\n",target);
                    //printf("evenIndex: %d\n",evenIndex);
                }
            }
            else if(index==-1)
            {
                
                //printf("算法进入了2:\n");
                if(oddIndex==0){oddIndextag=0;target=nneven_IOArray[0].startLpos-1;;}
                if(evenIndex!=0){xx=0;}
                else
                {
                    target=nnodd_IOArray[nonNegative(oddIndex-1)].startLpos+1;
                }
            }
             }
             
        }
        
    }

//#if DEBUG_LEVEL >= 1
    // printf("nn1算法得到的io序列:\n");
    // for (uint32_t i = 0; i < input->ioVec.len; ++i) {
    //      printf("%d ", output->sequence[i]);
    // }
    // printf("\n");
//#endif

    // 释放副本内存
    free(nnodd_IOArray);
    free(nneven_IOArray);
    return 0;
}

// 多轮线性SCAN-线性搜索
int32_t NearestNeighbor2(const InputParam *input, OutputParam *output){    
    // 创建 ioArray 的副本，用于排序
    //IOUintnn *nn_IOArray = (IOUintnn *)malloc(input->ioVec.len * sizeof(IOUintnn));
    int len1=0;
    int len0=0;
    for (uint32_t i = 0; i < input->ioVec.len; ++i) 
    {
        if(input->ioVec.ioArray[i].wrap%2==1)len1++;
        else len0++;
    }
    IOUintnn *nnodd_IOArray = (IOUintnn *)malloc(len1* sizeof(IOUintnn));;
    IOUintnn *nneven_IOArray  = (IOUintnn *)malloc(len0 * sizeof(IOUintnn));;

    uint32_t oddIndex = 0;
    uint32_t evenIndex = 0;
    
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        if(input->ioVec.ioArray[i].wrap%2==1)
        {
            nnodd_IOArray[oddIndex].endLpos = input->ioVec.ioArray[i].endLpos+5;
            nnodd_IOArray[oddIndex].startLpos = input->ioVec.ioArray[i].startLpos+5;
            nnodd_IOArray[oddIndex].wrap=input->ioVec.ioArray[i].wrap;
            nnodd_IOArray[oddIndex].id=input->ioVec.ioArray[i].id;
            nnodd_IOArray[oddIndex].visit=0;
            oddIndex++;
        }
        else
        {
            nneven_IOArray[evenIndex].endLpos = input->ioVec.ioArray[i].endLpos+5;
            nneven_IOArray[evenIndex].startLpos = input->ioVec.ioArray[i].startLpos+5;
            nneven_IOArray[evenIndex].wrap=input->ioVec.ioArray[i].wrap;
            nneven_IOArray[evenIndex].id=input->ioVec.ioArray[i].id;
            nneven_IOArray[evenIndex].visit=0;
            evenIndex++;
        }
    }
    
    // 使用 qsort 对副本进行排序
    qsort(nnodd_IOArray, oddIndex, sizeof(IOUintnn), nn_compare);
    qsort(nneven_IOArray, evenIndex, sizeof(IOUintnn), nn_compare);
    //input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status;
    //for (uint32_t i = 0; i < len0; ++i)
    //printf("偶数第%d个节点%dstarlops: %d\n",i,nneven_IOArray[i].id,nneven_IOArray[i].startLpos);
    //for (uint32_t i = 0; i < len1; ++i)
    //printf("奇数第%d个节点%dstarlops: %d\n",i,nnodd_IOArray[i].id,nnodd_IOArray[i].startLpos);
    uint32_t seqIndex = 0;
    int target_wrap=input->headInfo.wrap;
    int target=input->headInfo.lpos+5;
    //printf("循环2:\n");

    // int index;

    while(seqIndex<input->ioVec.len)
    {   //printf("循环: %d \n",len0);

        for (uint32_t i = 0; i < len0; ++i)
        {   //printf("循环 %d :\n",i);
            // printf("nneven_IOArray[i].visit: %d :\n",nneven_IOArray[i].visit);
            // printf("nneven_IOArray[i].startLpos: %d :\n",nneven_IOArray[i].startLpos);
            // printf("nnodd_IOArray[len1-1].startLpos: %d :\n",nnodd_IOArray[len1-1].startLpos);
            if((nneven_IOArray[i].visit==0)&&(nneven_IOArray[i].startLpos>target))
            {   //printf("进入了判断 :\n");
                int index=i+1;
                while(index<=len0-1&&nneven_IOArray[index].startLpos==nneven_IOArray[i].startLpos)
                {     
                    if((nneven_IOArray[index].visit==0)&&abs((int)nneven_IOArray[index].wrap-target_wrap)<abs((int)nneven_IOArray[i].wrap-target_wrap))
                    {
                        i=index;
                        //printf("有重复wrap情况:  \n");
                    }
                    index++;
                }
                nneven_IOArray[i].visit=1;
                target=nneven_IOArray[i].endLpos;
                output->sequence[seqIndex++]=nneven_IOArray[i].id;
                target_wrap=nneven_IOArray[i].wrap;
            }
            if(i==len0-1)
            {
                target=nnodd_IOArray[len1-1].startLpos+1;
            }
        }
        //printf("偶数wrap循环结束:\n");
        for (int i = len1-1; i >= 0; i--)
        {   //printf("循环 %d :\n",i);
            if(nnodd_IOArray[i].visit==0&&nnodd_IOArray[i].startLpos<target)
            {   
                int index=i-1;
                while((index>= 0)&&nnodd_IOArray[index].startLpos==nnodd_IOArray[i].startLpos)
                {     
                    if((nnodd_IOArray[index].visit==0)&&abs((int)nnodd_IOArray[index].wrap-target_wrap)<abs((int)nnodd_IOArray[i].wrap-target_wrap))
                    {
                        i=index;
                    }
                    index--;
                }
                nnodd_IOArray[i].visit=1;
                target=nnodd_IOArray[i].endLpos;
                output->sequence[seqIndex++]=nnodd_IOArray[i].id;
                target_wrap=nnodd_IOArray[i].wrap;
            }
            if(i==0)
            {
                target=nneven_IOArray[0].startLpos-1;
            }
        }
    }
    
#if DEBUG_LEVEL >= 1
    printf("nn1算法得到的io序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
         printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放副本内存
    free(nnodd_IOArray);
    free(nneven_IOArray);
    return 0;
}
// 贪心插入构造
int32_t GreedyInsert(const InputParam *input, OutputParam *output){
    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点
#ifndef USING_REAL_TIME_COST
   /* 生成邻接矩阵 */
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    for (uint32_t i = 0; i < len; ++i){
        adjMat[i] = (int *)malloc(sizeof(int) * len);
    }
    getAdjMat(input, len, adjMat);
#endif

    /* 基于贪心插入构造初始解， 注意巡回路径从1开始编号*/
    int *Next = (int *)calloc(len + 1, sizeof(int)); //记录当前城市的下一个城市
    // 从虚拟结点出发, 再到磁头结点
    Next[len - 1] = len;
    for (int i = 1; i <= len-2; ++i){ // 每次找一个城市
        int minDeltaCost = INF;
        int insertPrev = -1;
        // 评估磁头往后的k个插入位置
        int prev = len; // 从磁头开始
        while(prev != 0){
            int deltaCost;
            if(Next[prev] == 0) // 末尾位置插入
#ifndef USING_REAL_TIME_COST
                deltaCost = adjMat[prev-1][i-1];
#else
                deltaCost = getCost(input, len, prev-1, i-1);
#endif
            else // 中间位置插入
#ifndef USING_REAL_TIME_COST
                deltaCost = adjMat[prev-1][i-1] + adjMat[i-1][Next[prev]-1] - adjMat[prev-1][Next[prev]-1];
#else
                deltaCost = getCost(input, len, prev-1, i-1) + getCost(input, len, i-1, Next[prev]-1) - getCost(input, len, prev-1, Next[prev]-1);
#endif
            if (deltaCost < minDeltaCost){
                minDeltaCost = deltaCost;
                insertPrev = prev;
            }
            prev = Next[prev];
        }
        // 当前城市插入得到的最优位置
        if(Next[insertPrev] == 0){ // 末尾位置插入
            Next[insertPrev] = i;
        }
        else{ // 中间位置插入
            Next[i] = Next[insertPrev];
            Next[insertPrev] = i;
        }
    }
    // 路径链表拷贝至输出数组
    int curr = Next[len], idx = 0;
    while(curr != 0){
        output->sequence[idx++] = curr;
        curr = Next[curr];
    }

#if DEBUG_LEVEL >= 1
    printf("贪心插入算法得到的io序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
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
    free(Next);

    return 0;
}

// 最大匹配(MaxMatchingByHarryZHRsolver)-随机拼接
int MaxMatching_Random_qsort(const InputParam *input, OutputParam *output) {
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点

    // 初始化图
    Graph graph;
    initGraph(&graph, len);

    // 构建邻接表
    getAdjList1(input, &graph);

    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);

    // 调用最大匹配算法
    int* Next = solveMaxMatchingByHarryZHR(&matcher, &graph)->arr;

    int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点

    int idx = 0;
    int loops=0;
    for (int i = len - 1; i >= 0; --i) {
        if (!visited[i]) {
             loops++;
            int current = i;
           // printf("环路: %d", current);
            if (current != len - 1 && current != len - 2) {
                output->sequence[idx++] = current + 1;
            }
            while (Next[current] != -1 && !visited[Next[current]]) {
                visited[current] = 1;
                current = Next[current];
                if (current != len - 1 && current != len - 2) {
                    output->sequence[idx++] = current + 1;
                }
                // printf(" -> %d", current);
            }
            visited[current] = 1;
             //printf("\n");
        }
    }
    printf("环路个数: %d\n", loops);
// #if DEBUG_LEVEL >= 1
//     printf("最大匹配算法得到的 io 序列:\n");
//     for (uint32_t i = 0; i < input->ioVec.len; ++i) {
//         printf("%d ", output->sequence[i]);
//     }
//     printf("\n");
// #endif

    // 释放资源
    free(visited);
    freeMaxMatchingByHarryZHR(&matcher);
    
    // 释放图的邻接表内存
    for (int i = 0; i < graph.n; i++) {
        Edge *e = graph.lists[i].head;
        while (e) {
            Edge *tmp = e;
            e = e->next;
            free(tmp);
        }
    }
    free(graph.lists);

    return 0;
}

// 最大匹配(MaxMatchingByHarryZHRsolver)-随机拼接
int MaxMatching_Random_heap(const InputParam *input, OutputParam *output) {
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点

    // 初始化图
    Graph graph;
    initGraph(&graph, len);

    // 构建邻接表
    getAdjList2(input, &graph);
    printf("构建邻接表完成\n");
    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);
    printf("initMaxMatchingByHarryZHR完成\n");
    // 调用最大匹配算法
    int* Next = solveMaxMatchingByHarryZHR(&matcher, &graph)->arr;
  
    int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点

    int idx = 0;
    int loops=0;
    for (int i = len - 2; i >= 0; --i) {
        if (!visited[i]) {
             loops++;
            int current = i;
            //printf("环路: %d", current);
            if (current != len - 1 && current != len - 2) {
                output->sequence[idx++] = current + 1;
            }
            int firstloopnode=current;
            while (Next[current] != -1 && !visited[Next[current]]) {
                visited[current] = 1;
                current = Next[current];
                if (current != len - 1 && current != len - 2&&current!=firstloopnode) {
                    output->sequence[idx++] = current + 1;
                }
                 //printf(" -> %d", current);
            }
            visited[current] = 1;
             //printf("\n");
        }
    }
    printf("环路个数: %d\n", loops);

// #if DEBUG_LEVEL >= 1
    // printf("最大匹配算法得到的 io 序列:\n");
    // for (uint32_t i = 0; i < input->ioVec.len; ++i) {
    //     printf("%d ", output->sequence[i]);
    // }
    // printf("\n");
// #endif

    // 释放资源
    free(visited);
        printf("free(visited)\n");
    freeMaxMatchingByHarryZHR(&matcher);
    printf("freeMaxMatchingByHarryZHR\n");
    // 释放图的邻接表内存
    for (int i = 0; i < graph.n; i++) {
        Edge *e = graph.lists[i].head;
        while (e) {
            Edge *tmp = e;
            e = e->next;
            free(tmp);
        }
    }
    free(graph.lists);

    return 0;
}

// 最大匹配(MaxMatchingByHarryZHRsolver)-贪心拼接
int MaxMatching_Greedy(const InputParam *input, OutputParam *output) {
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
    // 初始化图
    Graph graph;
    initGraph(&graph, len);

    // 构建邻接表
    getAdjList2(input, &graph);
    printf("构建邻接表完成\n");
    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);
    printf("initMaxMatchingByHarryZHR完成\n");
    // 调用最大匹配算法
    int* Next = solveMaxMatchingByHarryZHR(&matcher, &graph)->arr;
    printf("最大匹配算法完成\n");
    int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点

    //printf("形成的环路：\n");
    int idx=0;
    int loops=0;
    int current=len-2;//i=len-1,确保从磁头开始构建序列
    while( current!=-1 && visited[current] != 1) {//i=len-1,确保从磁头开始构建序列
        visited[current] = 1;//标记访问
        loops++;
        int firstloopnode=current;
        if(current!=len-1&&current!=len-2){//若当前节点不是磁头节点或者虚拟节点，则放入结果序列中
            output->sequence[idx++]=current+1;
        }
        while (Next[current] != -1 && !visited[Next[current]]) {//继续基于当前节点向后搜索，直到把同一回路（即本回路）的节点全部放入结果序列为止
            visited[current] = 1;//标记访问
            current = Next[current];//传递下标索引
            if(current!=len-1&&current!=len-2&&current!=firstloopnode){//放入结果序列
                output->sequence[idx++]=current+1;
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
                int j=Next[i];//i的下一个节点j才是要判断是否连接的点
                 //printf("mincost=%d\n",mincost);
                if(mincost>getCost(input, len, current, j)-getCost(input, len, i, j)){
                    mincost=getCost(input, len, current, j)-getCost(input, len, i, j);
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

#if DEBUG_LEVEL >= 1
    printf("最大匹配算法得到的 io 序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放 matcher 和 adjMat 的内存
    free(visited);
    freeMaxMatchingByHarryZHR(&matcher);
    // 释放图的邻接表内存
    for (int i = 0; i < graph.n; i++) {
        Edge *e = graph.lists[i].head;
        while (e) {
            Edge *tmp = e;
            e = e->next;
            free(tmp);
        }
    }
    free(graph.lists);

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
     if(input->ioVec.len <= 7500){
        getbaseline(input, output);
        // 使用 LKH_NN 算法
        return LKH_NN(input, output);
        // 使用 LKH_MM 算法
        // return LKH_MM(input, output);
        //使用 LKH_nn1 算法
        // return LKH_nn1(input, output);
     }else{
        // return LKH_MM(input, output);
        // return MaxMatching_Greedy(input, output);
        // return MaxMatching_Greedy(input, output);
        // return NearestNeighbor(input, output);
        return NearestNeighbor2(input, output);
     }
    // // 使用 LKH_MM 算法
    // //return LKH_MM(input, output);

    // 使用 LKH_GI 算法
    // return LKH_GI(input, output); 

    // 使用 LKH_NN 算法
    //return LKH_NN(input, output);

    // 使用 LKH_nn1 算法
    //return LKH_nn1(input, output);

    
    // 使用最近邻贪心算法 
    //return NearestNeighbor(input, output);

    // 直接构造最近邻贪心算法 
    //return NearestNeighbor1(input, output);

    // 线性扫描
    //return NearestNeighbor2(input, output);

    // 使用贪心插入算法
    // return GreedyInsert(input, output);

    // 使用最大匹配-随机拼接算法
    //return MaxMatching_Random_qsort(input, output);
     //return MaxMatching_Random_heap(input, output);
    // 使用最大匹配-贪心拼接算法
    //return MaxMatching_Greedy(input, output);

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
