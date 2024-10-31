
#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "../public.h"

#ifdef __cplusplus
extern "C" {
#endif

int32_t AlgorithmRun(const InputParam *input, OutputParam *output);
double MyGetTime(); /* 返回实际时间：秒 */
int GetObjectValue(const HeadInfo *start, const HeadInfo *end); /* 返回目标函数值 */

#ifdef __cplusplus
}
#endif

#endif  // ALGORITHM_H