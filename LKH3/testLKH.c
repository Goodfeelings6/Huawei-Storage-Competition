#include "LKHInterface.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int main(){
    /* 生成邻接矩阵 */
    int len = 10;
    // int distance[5][5] = { // 对称实例
    //     {0, 10, 15, 20, 25},
    //     {10, 0, 35, 25, 30},
    //     {15, 35, 0, 30, 20},
    //     {20, 25, 30, 0, 15},
    //     {25, 30, 20, 15, 0}
    // };
    int distance[10][10] = { // 非对称实例
        {0, 12, 25, 30, 18, 40, 35, 20, 28, 15},
        {14, 0, 22, 29, 40, 25, 30, 21, 33, 19},
        {18, 30, 0, 20, 27, 45, 32, 40, 35, 22},
        {25, 28, 22, 0, 35, 20, 25, 38, 27, 30},
        {30, 35, 40, 18, 0, 15, 28, 22, 25, 37},
        {20, 40, 35, 32, 28, 0, 12, 24, 29, 18},
        {29, 34, 38, 27, 30, 20, 0, 17, 22, 25},
        {24, 32, 29, 20, 15, 18, 22, 0, 30, 26},
        {15, 28, 25, 18, 22, 32, 35, 30, 0, 14},
        {35, 40, 31, 29, 20, 25, 30, 28, 15, 0}
    };

    int **adjMat = (int **)malloc(sizeof(int *) * len);
    for (int i = 0; i < len; ++i){
        adjMat[i] = (int *)malloc(sizeof(int) * len);
    }
    assert(adjMat != NULL);
    for (int i = 0; i < len; ++i) {
        for (int j = 0; j < len; ++j) { 
            adjMat[i][j] = distance[i][j];
        }
    }

    /* 定义结果返回 */
    int *tourResult = (int *)malloc(sizeof(int) * len);
    assert(tourResult != NULL);
    int tourCost = 0;

    /* 确定LKH输入结构体 */
    LKHInput *lkhInput = (LKHInput *)malloc(sizeof(LKHInput));
    lkhInput->adjMat = adjMat;
    lkhInput->matDimension = len;
    lkhInput->intialTour = 0;
    lkhInput->fixEdge = 0;
    lkhInput->fixEdgeLen = 0; 
    lkhInput->lkhParameters = (LKHParameters *)malloc(sizeof(LKHParameters));
    loadDefaultParam(lkhInput->lkhParameters);
    /* 确定LKH输出结构体 */
    LKHOutput *lkhOutput = (LKHOutput *)malloc(sizeof(LKHOutput));
    lkhOutput->tourCost = 0;
    lkhOutput->tourResult = (int *)malloc(len * sizeof(int));

    /* 求解 */
    solveTSP(lkhInput, lkhOutput);


    printf("tourResult: ");
    for (int i = 0; i < len; ++i) {
        printf("%d ", lkhOutput->tourResult[i]);
    }
    printf("\ntourCost: %d\n", lkhOutput->tourCost);
    
    free(adjMat);
    free(tourResult);
    free(lkhInput->lkhParameters);
    free(lkhInput);
    free(lkhOutput);

    return 0;
}