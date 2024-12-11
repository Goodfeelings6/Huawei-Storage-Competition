#include "algorithm.h"
#include "MaxMatchingByHarryZHR.h"

// 最大匹配(MaxMatchingByHarryZHRsolver)-贪心拼接
// int MaxMatching_Greedy(const InputParam *input, OutputParam *output) {
//     uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
//     // 初始化图
//     Graph graph;
//     initGraph(&graph, len);

//     // 构建邻接表
//     getAdjList(input, &graph);

//     printf("构建邻接表完成\n");

//     // 释放图的邻接表内存
//     // for (int i = 0; i < graph.n; i++) {
//     //     Edge *e = graph.lists[i].head;
//     //     printf("node%d的邻接表；\n",i);
//     //     while (e) {
//     //         Edge *tmp = e;
//     //         printf("Node: %d, Cost: %d\n",e->target, e->cost);
//     //         e = e->next;
//     //         //free(tmp);
//     //     }
//     // }
//     // 初始化 MaxMatchingByHarryZHR 结构体
//     MaxMatchingByHarryZHR matcher;
//     initMaxMatchingByHarryZHR(&matcher, len);
//     printf("initMaxMatchingByHarryZHR完成\n");
//     // 调用最大匹配算法
//     int* Next = solveMaxMatchingByHarryZHR(&matcher, &graph)->arr;
//     printf("最大匹配算法完成\n");
//     int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点

//     // int* visited2 = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点

//     // int idx2 = 0;
//     // int loops2=0;
//     // for (int i = len - 2; i >= 0; --i) {
//     //     if (!visited2[i]) {
//     //          loops2++;
//     //         int current2 = i;
//     //         printf("环路: %d", current2);

//     //         int firstloopnode=current2;
//     //         while (Next[current2] != -1 && !visited2[Next[current2]]) {
//     //             visited[current2] = 1;
//     //             current2 = Next[current2];
               
//     //              printf(" -> %d", current2);
//     //         }
//     //         visited2[current2] = 1;
//     //          printf("\n");
//     //     }
//     // }
//     // printf("环路个数: %d\n", loops2);



//     printf("形成的环路：\n");
//     int idx=0;
//     int loops=0;
//     int current=len-2;//i=len-1,确保从磁头开始构建序列
//     while( current!=-1 && visited[current] != 1) {//i=len-1,确保从磁头开始构建序列
//         visited[current] = 1;//标记访问
//         loops++;
//         int firstloopnode=current;
//         if(current!=len-1&&current!=len-2){//若当前节点不是磁头节点或者虚拟节点，则放入结果序列中
//             output->sequence[idx++]=current+1;
//         }
//         while (Next[current] != -1 && !visited[Next[current]]) {//继续基于当前节点向后搜索，直到把同一回路（即本回路）的节点全部放入结果序列为止
//             visited[current] = 1;//标记访问
//             current = Next[current];//传递下标索引
//             if(current!=len-1&&current!=len-2&&current!=firstloopnode){//放入结果序列
//                 output->sequence[idx++]=current+1;
//             }
//         }
//         visited[current] = 1;
//         //放完本回路节点之后，遍历所有剩余节点，找到连接路径最短的其他回路节点(这里最短要考虑删边和加边)
//         int mincost=INF;
//         int nextnode=-1;//记录连接到的下一个回路的节点
//         int outtimes=0;
//         int intimes=0;
//         for (int i = len-3; i >= 0; --i) {
//             outtimes++;
//             if (!visited[i] ) {
//                 intimes++;
//                 int j=Next[i];//i的下一个节点j才是要判断是否连接的点
//                 if(mincost>getCost(input, len, current, j)-getCost(input, len, i, j)){
//                     mincost=getCost(input, len, current, j)-getCost(input, len, i, j);
//                     nextnode=j;
//                 }
//             }
//         }
//         current=nextnode;//把找到的连接路径最短的节点作为下一轮的起始节点
//     }
//     printf("环路个数: %d\n", loops);
    
// #if DEBUG_LEVEL >= 1
//     printf("最大匹配算法得到的 io 序列:\n");
//     for (uint32_t i = 0; i < input->ioVec.len; ++i) {
//         printf("%d ", output->sequence[i]);
//     }
//     printf("\n");
// #endif

//     // 释放 matcher 和 adjMat 的内存
//     free(visited);
//     freeMaxMatchingByHarryZHR(&matcher);
//     // 释放图的邻接表内存
//     for (int i = 0; i < graph.n; i++) {
//         Edge *e = graph.lists[i].head;
//         while (e) {
//             Edge *tmp = e;
//             e = e->next;
//             free(tmp);
//         }
//     }
//     free(graph.lists);

//     return 0;
// }


// 最大匹配(MaxMatchingByHarryZHRsolver)-贪心拼接
int MaxMatching_Greedy(const InputParam *input, OutputParam *output) {
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
// #ifndef USING_REAL_TIME_COST
   /* 生成邻接矩阵 */
    int **adjMat = (int **)malloc(sizeof(int *) * len);
    for (uint32_t i = 0; i < len; ++i){
        adjMat[i] = (int *)malloc(sizeof(int) * len);
    }
    getAdjMat(input, len, adjMat);

    // 矩阵转换
    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            adjMat[i][j] =  -adjMat[i][j];
        }
    }
// #endif

    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);
    int* Next = solveMaxMatchingByHarryZHR(&matcher, adjMat)->arr;
    int* visited = (int*)calloc(len, sizeof(int)); // 标记已经到达过的节点
    
    //printf("形成的环路：\n");
    int idx=0;
    int loops=0;
    int current=len-1;//i=len-1,确保从磁头开始构建序列
    while( current!=-1 && visited[current] != 1) {//i=len-1,确保从磁头开始构建序列
        visited[current] = 1;//标记访问
        loops++;
        if(current!=len-1&&current!=len-2){//若当前节点不是磁头节点或者虚拟节点，则放入结果序列中
            output->sequence[idx++]=current+1;
        }
        while (Next[current] != -1 && !visited[Next[current]]) {//继续基于当前节点向后搜索，直到把同一回路（即本回路）的节点全部放入结果序列为止
            visited[current] = 1;//标记访问
            current = Next[current];//传递下标索引
            if(current!=len-1&&current!=len-2){//放入结果序列
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
    for (uint32_t i = 0; i < len; ++i) {
        free(adjMat[i]);
    }
    free(adjMat);

    return 0;
}
