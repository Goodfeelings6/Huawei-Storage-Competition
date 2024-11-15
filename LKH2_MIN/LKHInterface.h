#ifndef _LKHInterface_H
#define _LKHInterface_H
#include "LKH.h"
#include "Genetic.h"

// LKH的参数接口
typedef struct{
    // 新增
    double TimeSpan; /* 统计改进值的时间跨度 */
    double PenaltyScoreInSecond; /* 20s后每多算1秒实际罚分 (相对Cost) */
    double ScheduleScoreInSecond; /* 20s内每少算1秒实际加分 (相对Cost) */
    const InputParam *OriginInput; /* 原始输入结构体 */
    
    char *ProblemFileName; /* PROBLEM_FILE */
    char *PiFileName; /* PI_FILE */
    char *InputTourFileName; /* INPUT_TOUR_FILE */
    char *OutputTourFileName;
    //char *SubproblemTourFileName;   /* SUBPROBLEM_TOUR_FILE */
    char *TourFileName; /* TOUR_FILE */
    int CandidateFiles; /* CANDIDATE_FILE */
    int MergeTourFiles; /* MERGE_TOUR_FILE */
    //int EdgeFiles;

    int AscentCandidates; /* ASCENT_CANDIDATES */
    int BackboneTrials;
    int Backtracking; /* BACKBONE_TRIALS */
    int CandidateSetSymmetric;
    int CandidateSetType; /* CANDIDATE_SET_TYPE */
    CrossoverFunction Crossover;
    int DelaunayPartitioning;
    int DelaunayPure;
    double Excess; /* EXCESS */
    int ExtraCandidates; /* EXTRA_CANDIDATES */
    int ExtraCandidateSetSymmetric; 
    int ExtraCandidateSetType; /* EXTRA_CANDIDATE_SET_TYPE */
    int Gain23Used; /* GAIN23 */
    int GainCriterionUsed; /* GAIN_CRITERION */
    double GridSize;
    int InitialPeriod; /* INITIAL_PERIOD */
    int InitialStepSize; /* INITIAL_STEP_SIZE */
    int InitialTourAlgorithm; /* INITIAL_TOUR_ALGORITHM */
    double InitialTourFraction; /* INITIAL_TOUR_FRACTION */
    int KarpPartitioning;
    int KCenterPartitioning;
    int KMeansPartitioning;
    int Kicks;  /* KICKS */
    int KickType; /* KICK_TYPE */
    int MaxBreadth; /* MAX_BREADTH */
    int MaxCandidates; /* MAX_CANDIDATES */
    int MaxPopulationSize;
    int MaxSwaps; /* MAX_SWAPS */
    int MaxTrials; /* MAX_TRIALS */
    int MoorePartitioning;
    int MoveType; /* MOVE_TYPE */
    int NonsequentialMoveType; /* NONSEQUENTIAL_MOVE_TYPE */
    GainType Optimum; /* OPTIMUM */
    int PatchingA; /* PATCHING_A */
    int PatchingC; /* PATCHING_C */
    int PatchingAExtended;
    int PatchingARestricted;
    int PatchingCExtended;
    int PatchingCRestricted;
    int Precision; /* PRECISION */
    int POPMUSIC_InitialTour; /* POPMUSIC_INITIAL_TOUR */
    int POPMUSIC_MaxNeighbors; /* POPMUSIC_MAX_NEIGHBORS */
    int POPMUSIC_SampleSize; /* POPMUSIC_SAMPLE_SIZE */
    int POPMUSIC_Solutions; /* POPMUSIC_SOLUTIONS */
    int POPMUSIC_Trials; /* POPMUSIC_TRIALS */
    int Recombination; /* RECOMBINATION */
    int RestrictedSearch; /* RESTRICTED_SEARCH */
    int RohePartitioning;
    int Runs; /* RUNS */
    unsigned int Seed; /* SEED */
    int SierpinskiPartitioning;
    int StopAtOptimum; /* STOP_AT_OPTIMUM */
    int Subgradient; /* SUBGRADIENT */
    int SubproblemBorders;
    int SubproblemsCompressed;
    int SubproblemSize;
    int SubsequentMoveType; /* SUBSEQUENT_MOVE_TYPE */
    int SubsequentPatching; /* SUBSEQUENT_PATCHING */
    double TimeLimit; /* TIME_LIMIT */
    double TotalTimeLimit; /* TOTAL_TIME_LIMIT */
    int TraceLevel; /* TRACE_LEVEL */
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
    double scheduleStartTime;    // 调度起始时间
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