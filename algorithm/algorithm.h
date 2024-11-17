
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "../public.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t AlgorithmRun(const InputParam *input, OutputParam *output);
double MyGetTime(); /* 返回实际时间：秒 */
int GetObjectValue(const HeadInfo *start, const HeadInfo *end); /* 返回目标函数值 */
int GetObjectValue_LKH(const HeadInfo *start, const HeadInfo *end); /* 返回目标函数值(LKH内部调用) */
extern double alpha;              // 读时延权重
extern double beta;               // 带体磨损权重
extern double gama;               // 电机磨损权重
extern double base_totaltime;     // 基线读时延
extern double base_tapeBeltWear;  // 基线带体磨损
extern double base_tapeMotorWear; // 基线电机磨损
extern double maxbase;            // 缩放系数

#ifdef __cplusplus
}
#endif

#endif  // ALGORITHM_H