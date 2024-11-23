
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <sys/time.h>
#include "../public.h"
#include "LKHInterface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ################## 宏定义 ################## */
#define USING_REAL_TIME_COST  // 是否启用实时的方式算代价，未定义则使用邻接矩阵方式

#define DEBUG_LEVEL 0
//#define DEBUG_LEVEL 1
// #define DEBUG_LEVEL 2

#define INF INT_MAX/200
#define HEAP_SIZE 50  // 堆的大小

/* ################## 全局变量声明 ################## */
extern double alpha;              // 读时延权重
extern double beta;               // 带体磨损权重
extern double gama;               // 电机磨损权重
extern double base_totaltime;     // 基线读时延
extern double base_tapeBeltWear;  // 基线带体磨损
extern double base_tapeMotorWear; // 基线电机磨损
extern double maxbase;            // 缩放系数（放大，减少精度丢失）
extern double final_alpha; // 最终权重
extern double final_beta;
extern double final_gama;

/* ################## 数据结构定义 ################## */
typedef struct
{
    uint32_t id;         // IO序号
    uint32_t wrap;       // 起始wrap
    uint32_t startLpos;  // 起始lpos
    uint32_t endLpos;    // 结束lpos
    uint32_t visit;      // 是否被选
} IOUintSS;

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
typedef struct {
    Edge *head;  // 邻接表的头指针
} AdjList;
typedef struct {
    AdjList *lists;  // 工人（左侧节点）对应的邻接表数组
    int n;  // 图中工人数量（也是任务数量）
} Graph;

/* ################## 函数声明 ################## */
/* algotithm.c */
// 设置需要调整的 LKH 的参数
void loadUserChangedParam(int matDimension, LKHParameters *p, double scheduleStartTime);
int32_t AlgorithmRun(const InputParam *input, OutputParam *output);
// 返回实际时间：秒
double MyGetTime();

/* GetCost.c */
// 返回目标函数值:两点间距离
int GetObjectValue(const HeadInfo *start, const HeadInfo *end); 
// 返回目标函数值:两点间距离(LKH内部调用)
int GetObjectValue_LKH(const HeadInfo *start, const HeadInfo *end); 
// 实时算代价
int getCost(const InputParam *input, int len, int i, int j);
// 邻接矩阵算代价
void getAdjMat(const InputParam *input, int len, int **adjMat);
// 邻接表算代价
void getAdjList(const InputParam *input, Graph *graph);

/* GetBaseline.c */
// 调用SCAN基线并分析IO序列以获取权重控制系数
void GetBaseline(const InputParam *input, OutputParam *output);

/* Sort.c */
// Sort算法
int32_t Sort(const InputParam *input, OutputParam *output);

/* Scan.c */
// Scan算法
int32_t Scan(const InputParam *input, OutputParam *output);

/* NearestNeighbor.c */
// 最近邻贪心构造
int32_t NearestNeighbor(const InputParam *input, OutputParam *output);

/* GreedyInsert.c */
// 贪心插入构造
int32_t GreedyInsert(const InputParam *input, OutputParam *output);

/* MaxMatching_Greedy.c */
// 最大匹配(MaxMatchingByHarryZHRsolver)+贪心拼接 构造
int MaxMatching_Greedy(const InputParam *input, OutputParam *output);

/* MultiSortScan.c */
// 多轮线性SCAN构造-线性搜索
int32_t MultiSortScan(const InputParam *input, OutputParam *output);

/* LKH_NN.c */
// LKH 算法 + 最近邻(NN)生成初始解
int32_t LKH_NN(const InputParam *input, OutputParam *output);

/* LKH_GI.c */
// LKH 算法 + 贪心插入(GI)生成初始解
int32_t LKH_GI(const InputParam *input, OutputParam *output);

/* LKH_MM.c */
// LKH 算法 + 最大匹配(MM)生成初始解
int32_t LKH_MM(const InputParam *input, OutputParam *output);

/* LKH_SS.c */
// LKH 算法 + 多轮线性SCAN(SS)生成初始解
int32_t LKH_SS(const InputParam *input, OutputParam *output);

#ifdef __cplusplus
}
#endif

#endif  // ALGORITHM_H