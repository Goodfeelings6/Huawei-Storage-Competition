#include "algorithm.h"

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
