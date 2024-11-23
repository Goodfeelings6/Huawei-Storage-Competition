#include "algorithm.h"

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
