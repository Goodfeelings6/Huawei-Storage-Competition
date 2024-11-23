#include "algorithm.h"
#include "MaxMatchingByHarryZHR.h"

// 最大匹配(MaxMatchingByHarryZHRsolver)-贪心拼接
int MaxMatching_Greedy(const InputParam *input, OutputParam *output) {
    uint32_t len = output->len + 2; // 增加了一个虚拟节点和磁头起始位置节点
    // 初始化图
    Graph graph;
    initGraph(&graph, len);

    // 构建邻接表
    getAdjList(input, &graph);
    // printf("构建邻接表完成\n");
    // 初始化 MaxMatchingByHarryZHR 结构体
    MaxMatchingByHarryZHR matcher;
    initMaxMatchingByHarryZHR(&matcher, len);
    // printf("initMaxMatchingByHarryZHR完成\n");
    // 调用最大匹配算法
    int* Next = solveMaxMatchingByHarryZHR(&matcher, &graph)->arr;
    // printf("最大匹配算法完成\n");
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

    // 释放 matcher 的内存
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
