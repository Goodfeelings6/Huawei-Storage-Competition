#ifndef SMART_SZX_GOAL_MAX_MATCHING_BY_HARRYZHR_H
#define SMART_SZX_GOAL_MAX_MATCHING_BY_HARRYZHR_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "Arr.h"  // 假设已经将 Arr.h 文件转换成了纯 C 的二维数组实现

#define INF 100000000

typedef int ID;
typedef int Cost;
typedef Cost (*CalcCost)(ID x, ID y, int **adjMat);  // 函数指针类型，用于计算成本矩阵中的值

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

int bfs(MaxMatchingByHarryZHR *matcher, int **adjMat) {
    while (!isLoopQueueEmpty(&matcher->q)) {
        ID u = frontLoopQueue(&matcher->q);
        popLoopQueue(&matcher->q);
        if (getArr2D(&matcher->visx, u, 0)) continue;
        setArr2D(&matcher->visx, u, 0, 1);
        for (ID i = 0; i < matcher->n; ++i) {
            if (adjMat[u][i] <= -INF) continue;
            if (getArr2D(&matcher->visy, i, 0)) continue;
            Cost delta = getArr2D(&matcher->lx, u, 0) + getArr2D(&matcher->ly, i, 0) - adjMat[u][i];
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

void KM(MaxMatchingByHarryZHR *matcher, int **adjMat) {

    for (ID i = 0; i < matcher->n; ++i) {

        setArr2D(&matcher->ly, i, 0, -INF);
        for (ID j = 0; j < matcher->n; ++j) {

            int value = adjMat[j][i];

            if (value > getArr2D(&matcher->ly, i, 0)) {
                setArr2D(&matcher->ly, i, 0, value);
              
            }
        }
        setArr2D(&matcher->lx, i, 0, 0);
    }
    for (ID i = 0; i < matcher->n; ++i) {
        clearMaxMatching(matcher);
        pushLoopQueue(&matcher->q, i);
        while (!bfs(matcher, adjMat));
    }
}

// 调用 KM 算法求解最大匹配问题
Arr2D* solveMaxMatchingByHarryZHR(MaxMatchingByHarryZHR *matcher, int **adjMat) {

    
    KM(matcher, adjMat);
    return &matcher->prx;
}

#endif // SMART_SZX_GOAL_MAX_MATCHING_BY_HARRYZHR_H
