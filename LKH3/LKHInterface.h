#ifndef _LKHInterface_H
#define _LKHInterface_H
#include "LKH.h"
#include "Genetic.h"
// LKH的参数接口
typedef struct{
    char *ProblemFileName;
    char *PiFileName;
    char *InputTourFileName;
    char *OutputTourFileName;
    char *TourFileName;
    int CandidateFiles;
    int MergeTourFiles;

    int AscentCandidates;
    int BackboneTrials;
    int Backtracking;
    int BWTSP_B;
    int BWTSP_Q;
    int BWTSP_L;
    int CandidateSetSymmetric;
    int CandidateSetType;
    CrossoverFunction Crossover;
    int DelaunayPartitioning;
    int DelaunayPure;
    int DemandDimension;
    double DistanceLimit;
    double Excess;
    int ExternalSalesmen;
    int ExtraCandidates;
    int ExtraCandidateSetSymmetric;
    int ExtraCandidateSetType;
    int Gain23Used;
    int GainCriterionUsed;
    double GridSize;
    int InitialPeriod;
    int InitialStepSize;
    int InitialTourAlgorithm;
    double InitialTourFraction;
    int KarpPartitioning;
    int KCenterPartitioning;
    int KMeansPartitioning;
    int Kicks;
    int KickType;
    int MaxBreadth;
    int MaxCandidates;
    int MaxPopulationSize;
    int MaxSwaps;
    int MaxTrials;
    int MoorePartitioning;
    int MoveType;
    int MoveTypeSpecial;
    int MTSPDepot;
    int MTSPMinSize;
    int MTSPMaxSize;
    int MTSPObjective;
    int NonsequentialMoveType;
    GainType Optimum ;
    int PatchingA;
    int PatchingC;
    int PatchingAExtended;
    int PatchingARestricted;
    int PatchingCExtended;
    int PatchingCRestricted;
    int Precision;
    int Probability;
    int POPMUSIC_InitialTour;
    int POPMUSIC_MaxNeighbors;
    int POPMUSIC_SampleSize;
    int POPMUSIC_Solutions;
    int POPMUSIC_Trials;
    int Recombination;
    int RestrictedSearch;
    int RohePartitioning;
    int Runs;
    int Salesmen;
    int Scale;
    unsigned int Seed;
    int SierpinskiPartitioning;
    int StopAtOptimum;
    int Subgradient;
    int SubproblemBorders;
    int SubproblemsCompressed;
    int SubproblemSize;
    int SubsequentMoveType;
    int SubsequentMoveTypeSpecial;
    int SubsequentPatching;
    double TimeLimit;
    double TotalTimeLimit;
    int TraceLevel;
    int TSPTW_Makespan;
    int MaxMatrixDimension;
} LKHParameters;
// 给参数结构体导入默认值
void loadDefaultParam(LKHParameters *p);

/* LKH的输入数据接口 */
typedef struct{
    int **adjMat;      // 邻接矩阵
    int matDimension;  // 矩阵维数，即结点数量
    int *intialTour;   // 初始解
    int *fixEdge;      // 固定边数组
    int fixEdgeLen;    // 固定边数量
    LKHParameters *lkhParameters; // LKH的参数
} LKHInput;

/* LKH的输出数据接口 */
typedef struct{
    int *tourResult;   // 结果路径存放数组
    long long tourCost;     // 结果路径的总成本
} LKHOutput;

/* Function prototypes: */
void ReadParameters(LKHInput* lkhInput); /* 参数设定 */
void ReadProblem(LKHInput* lkhInput); /* 问题信息设定、数据读入 */
void OutputTourResult(LKHOutput* lkhOutput); /* 提取求解结果 */
void ReSetLKH(); /* 重置 LKH内核 状态： 释放内存、全局变量重置、静态变量重置等 */
/**
 * @brief 求解TSP接口，即会调用 LKH内核 来求解
 * @param [in] lkhInput LKH输入
 * @param [out] lkhOutput LKH输出
 * @return int 
 */
int solveTSP(LKHInput* lkhInput, LKHOutput* lkhOutput);

#endif