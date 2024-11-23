#ifndef SMART_SZX_GOAL_MAX_MATCHING_ARR_H
#define SMART_SZX_GOAL_MAX_MATCHING_ARR_H
#include<stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    AllBits0 = 0,
    AllBits1 = -1,
    SafeMaxInt = 0x3F
} ResetOption;

// 定义一维数组结构体
typedef struct {
    int *arr;     // 指向数组数据的指针
    int length;   // 数组长度
} Arr;

// 定义二维数组结构体
typedef struct {
    int *arr;     // 指向数组数据的指针
    int rows;     // 行数
    int cols;     // 列数
    int totalSize; // 总元素数
} Arr2D;

// 一维数组相关函数

// 初始化一维数组
void initArr(Arr *a, int length) {
    a->arr = (int *)malloc(length * sizeof(int));
    a->length = (a->arr != NULL) ? length : 0;
}

// 释放一维数组内存
void clearArr(Arr *a) {
    if (a->arr != NULL) {
        free(a->arr);
        a->arr = NULL;
        a->length = 0;
    }
}

// 重置一维数组内容
void resetArr(Arr *a, ResetOption option) {
    if (a->arr != NULL) {
        memset(a->arr, option, a->length * sizeof(int));
    }
}

// 获取一维数组大小
int sizeArr(const Arr *a) {
    return a->length;
}

// 判断一维数组是否为空
int isEmptyArr(const Arr *a) {
    return (a->length == 0);
}

// 访问一维数组元素
int getArr(const Arr *a, int index) {
    return a->arr[index];
}

void setArr(Arr *a, int index, int value) {
    a->arr[index] = value;
}

// 二维数组相关函数

// 初始化二维数组
void initArr2D(Arr2D *a, int rows, int cols) {
    a->arr = (int *)malloc(rows * cols * sizeof(int));
    a->rows = (a->arr != NULL) ? rows : 0;
    a->cols = (a->arr != NULL) ? cols : 0;
    a->totalSize = rows * cols;
}

// 释放二维数组内存
void clearArr2D(Arr2D *a) {
    if (a->arr != NULL) {
        free(a->arr);
        a->arr = NULL;
        a->rows = a->cols = a->totalSize = 0;
    }
}

// 重置二维数组内容
void resetArr2D(Arr2D *a, ResetOption option) {
    if (a->arr != NULL) {
        memset(a->arr, option, a->totalSize * sizeof(int));
    }
}

// 获取二维数组行数
int rowsArr2D(const Arr2D *a) {
    return a->rows;
}

// 获取二维数组列数
int colsArr2D(const Arr2D *a) {
    return a->cols;
}

// 判断二维数组是否为空
int isEmptyArr2D(const Arr2D *a) {
    return (a->totalSize == 0);
}

// 获取二维数组元素
int getArr2D(const Arr2D *a, int row, int col) {
    return a->arr[row * a->cols + col];
}

// 设置二维数组元素
void setArr2D(Arr2D *a, int row, int col, int value) {
    a->arr[row * a->cols + col] = value;
}

#endif // SMART_SZX_GOAL_MAX_MATCHING_ARR_H
