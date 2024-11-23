#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "algorithm.h"

// 权重控制系数
double alpha = 0.3; // 读时延权重
double beta = 0.5; // 带体磨损权重
double gama = 0.2; // 电机磨损权重
double base_totaltime = 0;     // 基线读时延
double base_tapeBeltWear = 0;  // 基线带体磨损
double base_tapeMotorWear = 0; // 基线电机磨损
double maxbase = 0;            // 缩放系数（放大，减少精度丢失）
double final_alpha = 0; // 最终权重
double final_beta = 0;
double final_gama = 0;

// 设置需要调整的 LKH 的参数
void loadUserChangedParam(int matDimension, LKHParameters *p, double scheduleStartTime){
    if(matDimension <= 102){ // 10,50,100
        // 默认值
    }
    else if(matDimension <= 1002){ // 1000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag1000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 10;
        p->TimeSpan = 0.2;
    }
    else if(matDimension <= 2002){ // 2000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag2000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 5;
        p->TimeSpan = 0.2;
    }
    else if(matDimension <= 5002){ // 5000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag5000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 4;
        p->TimeSpan = 0.3;
    }
    else{ // 10000
        p->Subgradient = 0;
        p->CandidateSetType = POPMUSIC;
        p->POPMUSIC_InitialTour = 1;
        //@flag10000
		p->POPMUSIC_SampleSize = 20;
		p->POPMUSIC_Solutions = 3;
        p->TimeSpan = 0.5;
    }
    p->Runs = 1;
    p->TraceLevel = 1;
    p->TimeLimit = DBL_MAX; // 由总时间、一定时间跨度内的改进值共同控制退出即可
    p->TotalTimeLimit = 40; // 最大允许运行时间
    p->PenaltyScoreInSecond = alpha*1000/base_totaltime + ((2-(matDimension-2)/10000.0)/200.0)*maxbase; /* 20s后每多算1秒实际罚分 (相对Cost) */
    p->ScheduleScoreInSecond = alpha*1000/base_totaltime + (((matDimension-2)/10000.0)/200.0)*maxbase;/* 20s内每少算1秒实际加分 (相对Cost) */
    // printf("p->ScheduleScoreInSecond = %lf\n", p->ScheduleScoreInSecond);
    // printf("p->PenaltyScoreInSecond = %lf\n", p->PenaltyScoreInSecond);
    p->MoveType = 3;
}

/**
 * @brief  算法接口
 * @param  input            输入参数
 * @param  output           输出参数
 * @return int32_t          返回成功或者失败，RETURN_OK 或 RETURN_ERROR
 */
int32_t IOScheduleAlgorithm(const InputParam *input, OutputParam *output)
{    
    //  if(input->ioVec.len <= 7500){
    //     GetBaseline(input, output);
    //     // 使用 LKH_NN 算法
    //     return LKH_NN(input, output);
    //     // 使用 LKH_MM 算法
    //     // return LKH_MM(input, output);
    //  }else{
    //     // GetBaseline(input, output);
    //     // return LKH_MM(input, output);
    //     // return MaxMatching_Greedy(input, output);
    //     // return LKH_NN(input, output);
    //     // return NearestNeighbor(input, output);
    //        return MultiSortScan(input, output);
    //  }

    // GetBaseline(input, output);

    // 使用 LKH_MM 算法
    // return LKH_MM(input, output);

    // 使用 LKH_NN 算法
    // return LKH_NN(input, output);

    // 使用 LKH_GI 算法
    // return LKH_GI(input, output); 

    // 使用 LKH_SS 算法
    // return LKH_SS(input, output); 

    // 使用最近邻贪心算法 
    // return NearestNeighbor(input, output);

    // 多轮线性SCAN构造-线性搜索
    return MultiSortScan(input, output);

    // 使用贪心插入算法
    // return GreedyInsert(input, output);

    // 使用最大匹配-贪心拼接算法
    // return MaxMatching_Greedy(input, output);

    // 使用 Sort 算法
    // return Sort(input, output);

    // 使用 Scan 算法
    // return Scan(input, output);
}

/**
 * @brief  算法运行的主入口
 * @param  input            输入参数
 * @param  output           输出参数
 * @return uint32_t          返回成功或者失败，RETURN_OK 或 RETURN_ERROR
 */
int32_t AlgorithmRun(const InputParam *input, OutputParam *output)
{
    int32_t ret;

    ret = IOScheduleAlgorithm(input, output);

    return ret;
}
