#include "algorithm.h"

/* 返回目标函数值:两点间距离 */
int GetObjectValue(const HeadInfo *start, const HeadInfo *end){
    // double value = SeekTimeCalculate(start, end) * final_alpha + BeltWearTimes(start, end, NULL) * final_beta  + MotorWearTimes(start, end) * final_gama;
    double value = SeekTimeCalculate(start, end) * final_alpha + MotorWearTimes(start, end) * final_gama;
    
    //printf("value=%lf\n",value);
    // double value = BeltWearTimes(start, end, NULL) * final_beta + MotorWearTimes(start, end)*final_gama;
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

/* 返回目标函数值:两点间距离(LKH内部调用) */
int GetObjectValue_LKH(const HeadInfo *start, const HeadInfo *end){
    double value = SeekTimeCalculate(start, end) * final_alpha + BeltWearTimes(start, end, NULL) * final_beta + MotorWearTimes(start, end) * final_gama;
    return (int)value;
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

// 实时算代价
int getCost2(const InputParam *input, int len, int i, int j){
    if(i == j)
        return INF;
    else if(i==len-2){ // 磁头节点
        if(j==len-1){ // 磁头节点到虚拟节点的代价为极大值
            return INF;
        }
        else{ // 磁头节点到其他节点的代价为对应寻址时间
            HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
            HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
            return GetObjectValue(&start, &end);
        }
    }
    else if(i==len-1){ // 虚拟节点
        if(j==len-2){ // 虚拟节点到磁头节点的代价为0
            return 0;
        }
        else{ // 虚拟节点到其他节点的代价为极大值
            return INF;
        }
    }
    else{ // 其他节点
        if(j==len-2){// 其他节点到磁头节点的代价为极大值
            return INF;
        }
        else if(j==len-1){// 其他节点到虚拟节点的代价为0
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

// // 邻接表算代价
// void getAdjList(const InputParam* input, Graph* graph) {
//     uint32_t len = graph->n;

//     for (uint32_t i = 0; i < len - 1; ++i) {
//         NodeCost heap[HEAP_SIZE];  // 大顶堆数组
//         int heapSize = 0;
//         int iwrap=input->ioVec.ioArray[i].wrap;
//         int index_samewrap_close=0;
//         int mindis=INF;
//         for (uint32_t j = 0; j < len - 1; ++j) {
//             int jwrap=input->ioVec.ioArray[j].wrap;
//             if (i == j) {
//                 addEdge(graph, i, i, -INF+1);
//                 continue;
//             }
          
//             if(iwrap==jwrap){//对同wrap的处理
//                 if(iwrap%2==0){//正向
//                     if(input->ioVec.ioArray[j].startLpos < input->ioVec.ioArray[i].endLpos){
//                         continue;
//                     }
//                     else if(mindis >= input->ioVec.ioArray[j].startLpos - input->ioVec.ioArray[i].endLpos){
//                         mindis=input->ioVec.ioArray[j].startLpos - input->ioVec.ioArray[i].endLpos;
//                         index_samewrap_close=j;//同wrap里i最近的IO块是j
//                     }
                    
//                 }else{
//                     if(input->ioVec.ioArray[j].startLpos > input->ioVec.ioArray[i].endLpos){
//                         continue;
//                     }
//                     else if(mindis >= input->ioVec.ioArray[i].endLpos -input->ioVec.ioArray[j].startLpos){
//                         mindis=input->ioVec.ioArray[i].endLpos -input->ioVec.ioArray[j].startLpos;
//                         index_samewrap_close=j;//同wrap里i最近的IO块是j
//                     }
//                 }
//             continue;
//             }
//             HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
//             HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
//             int cost = GetObjectValue(&start, &end);
//             //printf("%d->%d,cost:%d\n",i,j,cost);
//             NodeCost newElement = {j, cost};
//             insertMaxHeap(heap, &heapSize, newElement);
//         }

//         // 将堆中元素添加到图的邻接表中
//         //printf("node%d的邻接表；\n",i);
//         for (int k = 0; k < heapSize; ++k) {
//             addEdge(graph, i, heap[k].node, -heap[k].cost);
//            // printf("Node: %d, Cost: %d\n", heap[k].node, heap[k].cost);
        
//             // HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
//             // HeadInfo end = {input->ioVec.ioArray[index_samewrap_close].wrap, input->ioVec.ioArray[index_samewrap_close].startLpos, HEAD_RW};
//             // int cost = GetObjectValue(&start, &end);
//             // //if(cost<100) printf("cost:%d\n",cost);
//             // addEdge(graph, i, index_samewrap_close, -cost);
//         }
//     }

//     // 处理虚拟节点和磁头节点的代价
//     for (uint32_t i = 0; i < len - 1; ++i) {
//         addEdge(graph, i, len - 1, 0);

//         // HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
//         // HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
//         // int headCost = GetObjectValue(&start, &end);
//         // addEdge(graph, len - 1, i, -headCost);
//     }
//     addEdge(graph, len - 1, len - 2, 0);
// }

