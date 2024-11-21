#ifndef _PRIORITYQUEUE_H
#define _PRIORITYQUEUE_H

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MaxMatchingByHarryZHR.h" 

int PriorityQueuesize=50;

//定义邻接表结构体
typedef struct EdgeGC {
    int id;  // 右侧节点（任务）
    struct EdgeGC *next;  // 下一条边
} EdgeGC;

typedef struct {
    EdgeGC *head;  // 邻接表的头指针
} AdjListGC;

typedef struct {
    AdjListGC *lists;  // 工人（左侧节点）对应的邻接表数组
    int n;  // 图中工人数量（也是任务数量）
} GraphGC;

void initGraphGC(GraphGC *graph, int n) {
    graph->n = n;
    graph->lists = (AdjListGC*)malloc(n * sizeof(AdjListGC));
    for (int i = 0; i < n; i++) {
        graph->lists[i].head = NULL;  // 初始化每个工人的邻接表为空
    }
}

void addEdgeGC(GraphGC *graph, int u, int v) {
    EdgeGC *newEdge = (EdgeGC*)malloc(sizeof(EdgeGC));
    newEdge->id= v;  // 设置目标节点
    newEdge->next = graph->lists[u].head;  // 插入到链表头部
    graph->lists[u].head = newEdge; // 更新节点u的邻接表
}


// 插入到最小优先队列的函数
void insertMinPriorityQueue(NodeCost* queue, int* queueSize, NodeCost newElement) {
    // 如果队列未满，直接插入并排序
    if (*queueSize < PriorityQueuesize) {
        queue[*queueSize] = newElement; // 插入到末尾
        (*queueSize)++;

        // 插入后维护最小优先队列的顺序
        for (int i = *queueSize - 1; i > 0 && queue[i].cost < queue[i - 1].cost; i--) {
            // 交换顺序
            NodeCost temp = queue[i];
            queue[i] = queue[i - 1];
            queue[i - 1] = temp;
        }
    } else {
        // 队列已满，检查新元素是否小于队列中最大元素
        if (newElement.cost < queue[*queueSize - 1].cost) {
            // 替换掉最大元素
            queue[*queueSize - 1] = newElement;

            // 从尾部向前维护顺序
            for (int i = *queueSize - 1; i > 0 && queue[i].cost < queue[i - 1].cost; i--) {
                // 交换顺序
                NodeCost temp = queue[i];
                queue[i] = queue[i - 1];
                queue[i - 1] = temp;
            }
        }
    }
}
#endif