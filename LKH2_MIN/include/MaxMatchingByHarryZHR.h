#ifndef SMART_SZX_GOAL_MAX_MATCHING_BY_HARRYZHR_H
#define SMART_SZX_GOAL_MAX_MATCHING_BY_HARRYZHR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "Arr.h"  // 假设已经将 Arr.h 文件转换成了纯 C 的二维数组实现

#define INF INT_MAX/200
#define HEAP_SIZE 50  // 堆的大小
typedef int ID;
typedef double Cost;

// 结构体用于存储节点和代价
typedef struct {
    int node;
    int cost;
} NodeCost;

//定义邻接表结构体
typedef struct Edge {
    ID target;  // 右侧节点（任务）
    Cost cost;  // 边的权重（匹配成本）
    struct Edge *next;  // 下一条边
} Edge;

// 比较函数，用于排序
int compareEdges(const void *a, const void *b) {
    return ((NodeCost *)b)->cost - ((NodeCost *)a)->cost;
}

typedef struct {
    Edge *head;  // 邻接表的头指针
} AdjList;

typedef struct {
    AdjList *lists;  // 工人（左侧节点）对应的邻接表数组
    int n;  // 图中工人数量（也是任务数量）
} Graph;

void initGraph(Graph *graph, int n) {
    graph->n = n;
    graph->lists = (AdjList*)malloc(n * sizeof(AdjList));
    for (int i = 0; i < n; i++) {
        graph->lists[i].head = NULL;  // 初始化每个工人的邻接表为空
    }
}

void addEdge(Graph *graph, ID u, ID v, Cost cost) {
    Edge *newEdge = (Edge*)malloc(sizeof(Edge));
    newEdge->target = v;  // 设置目标节点
    newEdge->cost = cost;  // 设置代价
    newEdge->next = graph->lists[u].head;  // 插入到链表头部
    graph->lists[u].head = newEdge; // 更新节点u的邻接表
}


// 插入到最大堆中，堆满时如果新元素比堆顶小则替换堆顶元素
void insertMaxHeap(NodeCost* heap, int* heapSize, NodeCost newElement) {
    if (*heapSize < HEAP_SIZE) {
        // 堆未满，直接插入新元素并上浮调整
        heap[(*heapSize)++] = newElement;
        // 上浮调整堆
        for (int i = *heapSize - 1; i > 0 && heap[i].cost > heap[(i - 1) / 2].cost; i = (i - 1) / 2) {
            NodeCost temp = heap[i];
            heap[i] = heap[(i - 1) / 2];
            heap[(i - 1) / 2] = temp;
        }
    } else if (newElement.cost < heap[0].cost) {
        // 堆已满，且新元素比堆顶小，则替换堆顶
        heap[0] = newElement;
        // 下沉调整堆
        int i = 0;
        while (2 * i + 1 < HEAP_SIZE) {
            int largest = i;
            int left = 2 * i + 1;
            int right = 2 * i + 2;
            if (left < HEAP_SIZE && heap[left].cost > heap[largest].cost)
                largest = left;
            if (right < HEAP_SIZE && heap[right].cost > heap[largest].cost)
                largest = right;
            if (largest == i)
                break;
            NodeCost temp = heap[i];
            heap[i] = heap[largest];
            heap[largest] = temp;
            i = largest;
        }
    }
}

// // 打印堆的内容
// void printHeap(NodeCost* heap, int heapSize, int iteration, int node) {
//     printf("Heap after iteration %d for node %d:\n", iteration, node);
//     for (int i = 0; i < heapSize; ++i) {
//         printf("Node: %d, Cost: %d\n", heap[i].node, heap[i].cost);
//     }
//     printf("\n");
//}
// 定义循环队列 LoopQueue 结构体
typedef struct {
    ID *q;            // 队列数组
    int len;          // 队列大小
    int head;         // 队首索引
    int tail;         // 队尾索引
} LoopQueue;

// 初始化循环队列
void initLoopQueue(LoopQueue *queue, int size) {
    queue->q = (ID*)malloc(size * sizeof(ID));
    queue->len = size;
    queue->head = 0;
    queue->tail = 0;
}

// 清空循环队列
void clearLoopQueue(LoopQueue *queue) {
    queue->head = queue->tail = 0;
}

// 判断队列是否为空
int isLoopQueueEmpty(LoopQueue *queue) {
    return (queue->head == queue->tail);
}

// 获取队首元素
ID frontLoopQueue(LoopQueue *queue) {
    return queue->q[queue->head];
}

// 向队列添加元素
void pushLoopQueue(LoopQueue *queue, ID item) {
    queue->q[queue->tail] = item;
    queue->tail = (queue->tail + 1) % queue->len;
}

// 从队列弹出元素
void popLoopQueue(LoopQueue *queue) {
    queue->head = (queue->head + 1) % queue->len;
}

// 释放循环队列内存
void freeLoopQueue(LoopQueue *queue) {
    free(queue->q);
}

// MaxMatchingByHarryZHR 结构体
typedef struct {
    //CalcCost cost;      // 成本矩阵的函数指针
    int n;              // n 表示工人和任务的数量
    Arr2D lx;           // 左侧顶标二维数组
    Arr2D ly;           // 右侧顶标二维数组
    Arr2D slack;        // 松弛二维数组
    Arr2D prx;          // 左侧匹配结果二维数组
    Arr2D pry;          // 右侧匹配结果二维数组
    Arr2D pre;          // 前驱节点二维数组
    Arr2D visx;         // 左侧访问标记二维数组
    Arr2D visy;         // 右侧访问标记二维数组
    LoopQueue q;        // 循环队列
} MaxMatchingByHarryZHR;

// 初始化 MaxMatchingByHarryZHR 结构体
void initMaxMatchingByHarryZHR(MaxMatchingByHarryZHR *matcher, ID dimension) {

    matcher->n = dimension;
    
    // 初始化各个二维数组
    initArr2D(&matcher->lx, dimension, 1);
    initArr2D(&matcher->ly, dimension, 1);
    initArr2D(&matcher->slack, dimension, 1);
    initArr2D(&matcher->prx, dimension, 1);
    initArr2D(&matcher->pry, dimension, 1);
    initArr2D(&matcher->pre, dimension, 1);
    initArr2D(&matcher->visx, dimension, 1);
    initArr2D(&matcher->visy, dimension, 1);

    resetArr2D(&matcher->lx, AllBits1);
    resetArr2D(&matcher->ly, AllBits1);
    resetArr2D(&matcher->prx, AllBits1);
    resetArr2D(&matcher->pry, AllBits1);
    resetArr2D(&matcher->pre, AllBits1);

    initLoopQueue(&matcher->q, dimension);
}

// 释放 MaxMatchingByHarryZHR 内存
void freeMaxMatchingByHarryZHR(MaxMatchingByHarryZHR *matcher) {
    clearArr2D(&matcher->lx);
    clearArr2D(&matcher->ly);
    clearArr2D(&matcher->slack);
    clearArr2D(&matcher->prx);
    clearArr2D(&matcher->pry);
    clearArr2D(&matcher->pre);
    clearArr2D(&matcher->visx);
    clearArr2D(&matcher->visy);
    freeLoopQueue(&matcher->q);
}

// KM 算法核心实现
int check(MaxMatchingByHarryZHR *matcher, ID x) {
    setArr2D(&matcher->visy, x, 0, 1);
    if (getArr2D(&matcher->pry, x, 0) >= 0) {
        pushLoopQueue(&matcher->q, getArr2D(&matcher->pry, x, 0));
        return 0;
    }
    while (x >= 0) {
        setArr2D(&matcher->pry, x, 0, getArr2D(&matcher->pre, x, 0));
        ID temp = x;
        x = getArr2D(&matcher->prx, getArr2D(&matcher->pre, temp, 0), 0);
        setArr2D(&matcher->prx, getArr2D(&matcher->pre, temp, 0), 0, temp);
    }
    return 1;
}

void clearMaxMatching(MaxMatchingByHarryZHR *matcher) {
    resetArr2D(&matcher->visx, AllBits0);
    resetArr2D(&matcher->visy, AllBits0);
    clearLoopQueue(&matcher->q);
    resetArr2D(&matcher->slack, SafeMaxInt);
}

int bfs(MaxMatchingByHarryZHR *matcher, Graph *graph) {
    while (!isLoopQueueEmpty(&matcher->q)) {
        ID u = frontLoopQueue(&matcher->q);
        popLoopQueue(&matcher->q);
        if (getArr2D(&matcher->visx, u, 0)) continue;
        setArr2D(&matcher->visx, u, 0, 1);

        // 遍历工人 u 的邻接表
        for (Edge *e = graph->lists[u].head; e != NULL; e = e->next) {
            ID i = e->target;  // 任务 i
            Cost edgeCost = e->cost;  // 匹配成本
           if (edgeCost <= -INF) continue;
            if (getArr2D(&matcher->visy, i, 0)) continue;
            Cost delta = getArr2D(&matcher->lx, u, 0) + getArr2D(&matcher->ly, i, 0) - edgeCost;
            if (delta < getArr2D(&matcher->slack, i, 0)) {
                setArr2D(&matcher->slack, i, 0, delta);
                setArr2D(&matcher->pre, i, 0, u);
                if (!delta && check(matcher, i)) return 1;
            }
        }
    }
    Cost delta = INF;
    for (ID i = 0; i < matcher->n; ++i) {
        if (!getArr2D(&matcher->visy, i, 0)) {
            delta = (delta < getArr2D(&matcher->slack, i, 0)) ? delta : getArr2D(&matcher->slack, i, 0);
        }
    }
    for (ID i = 0; i < matcher->n; ++i) {
        if (getArr2D(&matcher->visx, i, 0)) setArr2D(&matcher->lx, i, 0, getArr2D(&matcher->lx, i, 0) - delta);
        if (getArr2D(&matcher->visy, i, 0)) setArr2D(&matcher->ly, i, 0, getArr2D(&matcher->ly, i, 0) + delta);
        else setArr2D(&matcher->slack, i, 0, getArr2D(&matcher->slack, i, 0) - delta);
    }
    for (ID i = 0; i < matcher->n; ++i) {
        if (!getArr2D(&matcher->visy, i, 0) && !getArr2D(&matcher->slack, i, 0) && check(matcher, i)) return 1;
    }

    return 0;
}

// int bfs(MaxMatchingByHarryZHR *matcher, Graph *graph) {
//     // 记录是否找到增广路径的标志
//     int foundAugmentingPath = 0;

//     while (!isLoopQueueEmpty(&matcher->q)) {
//         ID u = frontLoopQueue(&matcher->q);
//         popLoopQueue(&matcher->q);
        
//         if (getArr2D(&matcher->visx, u, 0)) continue;  // 已经访问过的工人跳过
//         setArr2D(&matcher->visx, u, 0, 1);  // 标记工人 u 已经访问

//         // 遍历工人 u 的邻接表
//         for (Edge *e = graph->lists[u].head; e != NULL; e = e->next) {
//             ID i = e->target;  // 任务 i
//             Cost edgeCost = e->cost;  // 匹配成本

//             if (getArr2D(&matcher->visy, i, 0)) continue;  // 任务 i 已经访问过，跳过

//             Cost delta = getArr2D(&matcher->lx, u, 0) + getArr2D(&matcher->ly, i, 0) - edgeCost;
//             if (delta < getArr2D(&matcher->slack, i, 0)) {
//                 setArr2D(&matcher->slack, i, 0, delta);
//                 setArr2D(&matcher->pre, i, 0, u);
                
//                 // 如果 delta == 0，尝试增广路径
//                 if (!delta && check(matcher, i)) {
//                     foundAugmentingPath = 1;  // 找到了增广路径
//                     break;  // 退出当前循环，结束 bfs
//                 }
//             }
//         }

//         if (foundAugmentingPath) {
//             break;  // 如果找到了增广路径，直接退出 bfs
//         }
//     }

//     // 如果没有增广路径，进行松弛操作
//     if (!foundAugmentingPath) {
//         Cost delta = INF;
//         for (ID i = 0; i < matcher->n; ++i) {
//             if (!getArr2D(&matcher->visy, i, 0)) {
//                 delta = (delta < getArr2D(&matcher->slack, i, 0)) ? delta : getArr2D(&matcher->slack, i, 0);
//             }
//         }

//         // 更新 lx 和 ly，松弛 slack 值
//         for (ID i = 0; i < matcher->n; ++i) {
//             if (getArr2D(&matcher->visx, i, 0)) setArr2D(&matcher->lx, i, 0, getArr2D(&matcher->lx, i, 0) - delta);
//             if (getArr2D(&matcher->visy, i, 0)) setArr2D(&matcher->ly, i, 0, getArr2D(&matcher->ly, i, 0) + delta);
//             else setArr2D(&matcher->slack, i, 0, getArr2D(&matcher->slack, i, 0) - delta);
//         }

//         // 如果没有增广路径，再次尝试增广路径
//         for (ID i = 0; i < matcher->n; ++i) {
//             if (!getArr2D(&matcher->visy, i, 0) && !getArr2D(&matcher->slack, i, 0) && check(matcher, i)) return 1;
//         }
//     }

//     return 0;
// }

void KM(MaxMatchingByHarryZHR *matcher, Graph *graph) {
    for (ID i = 0; i < matcher->n; ++i) {
        setArr2D(&matcher->ly, i, 0, -INF);
        
        // 遍历工人 i 的邻接表来更新 lx 和 ly
        for (Edge *e = graph->lists[i].head; e != NULL; e = e->next) {
            ID j = e->target;
            int value = e->cost;  // 匹配成本

            if (value > getArr2D(&matcher->ly, j, 0)) {
                setArr2D(&matcher->ly, j, 0, value);
            }
        }
        setArr2D(&matcher->lx, i, 0, 0);
    }
    for (ID i = 0; i < matcher->n; ++i) {    
        clearMaxMatching(matcher);
        pushLoopQueue(&matcher->q, i);
        while (!bfs(matcher, graph));
    }
}


// 调用 KM 算法求解最大匹配问题
// Arr2D* solveMaxMatchingByHarryZHR(MaxMatchingByHarryZHR *matcher, int **adjMat) {

    
//     KM(matcher, adjMat);
//     return &matcher->prx;
// }
Arr2D* solveMaxMatchingByHarryZHR(MaxMatchingByHarryZHR *matcher, Graph *graph) {
    KM(matcher, graph);
    return &matcher->prx;
}

#endif // SMART_SZX_GOAL_MAX_MATCHING_BY_HARRYZHR_H
