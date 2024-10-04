#ifndef _LKHInterface_H
#define _LKHInterface_H

/* Function prototypes: */

void ReadParameters(void); /* 参数设定 */
void ReadProblem(int **adjMat, int MatDimension); /* 问题信息设定、数据读入 */
void OutputTourResult(int* tourResult, int *tourCost); /* 提取求解结果 */
void ReSetLKH(); /* 重置 LKH内核 状态： 释放内存、全局变量重置、静态变量重置等 */
/**
 * @brief 求解TSP接口，即会调用 LKH内核(LKH2 算法)来求解
 * @param [in] adjMat 邻接矩阵
 * @param [in] MatDimension 矩阵维数，即结点数量
 * @param [out] tourResult 结果路径存放数组
 * @param [out] tourCost 结果路径的总成本
 * @return int 
 */
int solveTSP(int **adjMat, int MatDimension, int *tourResult, int *tourCost);

#endif