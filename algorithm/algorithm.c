#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "algorithm.h"
#include "LKHInterface.h"

/**
 * @brief  算法接口
 * @param  input            输入参数
 * @param  output           输出参数
 * @return int32_t          返回成功或者失败，RETURN_OK 或 RETURN_ERROR
 */
int32_t IOScheduleAlgorithm(const InputParam *input, OutputParam *output)
{    
    // 使用 LKH 算法
    // return LKH(input, output);
    // 使用 Sort 算法
    // return Sort(input, output);
    // 使用 Scan 算法
    return Scan(input, output);
}

// LKH 算法
int32_t LKH(const InputParam *input, OutputParam *output)
{    
    int32_t ret = 0;

    /* 生成邻接矩阵 */
    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点
    uint32_t **adjMat = (uint32_t **)malloc(sizeof(uint32_t *) * len);
    for (uint32_t i = 0; i < len; ++i){
        adjMat[i] = (uint32_t *)malloc(sizeof(uint32_t) * len);
    }
    assert(adjMat != NULL);
   
    for (uint32_t i = 0; i < len - 2; ++i) {//为了节点能和id号对应，还是将它们从0开始对应行列
        for (uint32_t j = 0; j < len - 2; ++j) {
            if(i == j)
                adjMat[i][j] = 0;
            else {
                HeadInfo start = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].endLpos, HEAD_RW};
                HeadInfo end = {input->ioVec.ioArray[j].wrap, input->ioVec.ioArray[j].startLpos, HEAD_RW};
                adjMat[i][j] = SeekTimeCalculate(&start, &end);
            }
        }
    }

     /*len-2行len-2列是虚拟节点，len-1行len-1列是磁头节点，为虚拟节点和磁头节点设置代价 */
    for (uint32_t i = 0; i < len - 2; ++i) {
        adjMat[len-2][i] = 0xfffff;  // 虚拟节点到其他节点的代价为极大值，则不可能达到其他点
        adjMat[i][len-2] = 0;  // 其他节点到虚拟节点的代价为0
        
        HeadInfo start = {input->headInfo.wrap, input->headInfo.lpos, input->headInfo.status};
        HeadInfo end = {input->ioVec.ioArray[i].wrap, input->ioVec.ioArray[i].startLpos, HEAD_RW};
        adjMat[len-1][i] = SeekTimeCalculate(&start, &end);  // 磁头节点到其他节点的代价为对应寻址时间
        adjMat[i][len-1] = 0xfffff;  // 其他节点到磁头节点的代价为极大值，则不可能返回该点
    }
    adjMat[len-2][len-1] = 0; // 虚拟节点到磁头节点的代价为0，到其他结点代价无穷，相当于固定了磁头节点为第二个节点
    adjMat[len-2][len-2] = 0; // 虚拟节点到自身的代价为0
    adjMat[len-1][len-1] = 0; // 磁头节点到自身的代价为0
    adjMat[len-1][len-2] = 0xfffff;  // 磁头节点到虚拟节点的代价为极大值

    /* 打印邻接矩阵 */
    for (uint32_t i = 0; i < len; ++i) {
        for (uint32_t j = 0; j < len; ++j) {
            if (adjMat[i][j] == 0xfffff) {
                printf("INF\t");  // 打印无穷大（不可达）的情况
            } else {
                printf("%u\t ", adjMat[i][j]);  // 打印代价
            }
        }
        printf("\n");
    }

    /* 求解 */
    int Cost = 0;
    int *result_seqence = (int *)malloc(len * sizeof(int));
    ret = solveTSP((int **)adjMat, len, result_seqence, &Cost);
   
    // 处理解
    int i, j;
    for (i = 0; i < len; ++i){
        if(result_seqence[i] == len) // i指向磁头结点
            break;
    }
    for (j = 0; j < len - 2; ++j){ // 排序结果赋值到输出
        if(i + 1 > len - 1)
            i = -1;
        output->sequence[j] = result_seqence[++i];
    }

    /* 打印求解结果：输出output->sequence的内容 */
    printf("len:%d\n",len);
    printf("Output sequence:\n");
    for (uint32_t i = 0; i < len - 2; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");

    /* 算法示例：先入先出算法 */
    // output->len = input->ioVec.len;
    // for (uint32_t i = 0; i < output->len; i++) {
    //     output->sequence[i] = input->ioVec.ioArray[i].id;
    // }

    /* 调用公共函数示例：调用电机寻址、带体磨损、电机磨损函数 */
    // HeadInfo start = {input->ioVec.ioArray[0].wrap, input->ioVec.ioArray[0].endLpos, HEAD_RW};
    // HeadInfo end = {input->ioVec.ioArray[1].wrap, input->ioVec.ioArray[1].endLpos, HEAD_RW};
    // int32_t seekT = 0;
    // int32_t beltW = 0;
    // int32_t motorW = 0;
    //    for (uint32_t i = 0; i < 10000; i++) {
    //        seekT = SeekTimeCalculate(&start, &end);
    //        beltW = BeltWearTimes(&start, &end, NULL);
    //        motorW = MotorWearTimes(&start, &end);
    //    }

    // /* 调用公共函数示例：调用IO读写时间函数 */
    // uint32_t rwT = ReadTimeCalculate(abs(input->ioVec.ioArray[0].endLpos - input->ioVec.ioArray[0].startLpos));

    return ret;
}

// Sort 算法的自定义比较函数：首先按 wrap 排序，相同 wrap 的按 startLpos 排序
int Sort_compare(const void *a, const void *b) {
    IOUint *ioA = (IOUint *)a;
    IOUint *ioB = (IOUint *)b;

    if (ioA->wrap != ioB->wrap) {
        return (ioA->wrap > ioB->wrap) - (ioA->wrap < ioB->wrap);  // 按 wrap 升序
    } else {
        return (ioA->startLpos > ioB->startLpos) - (ioA->startLpos < ioB->startLpos);  // 按 startLpos 升序
    }
}

int32_t Sort(const InputParam *inputParam, OutputParam *output){
    // 创建 ioArray 的副本，用于排序
    IOUint *sort_IOArray = (IOUint *)malloc(inputParam->ioVec.len * sizeof(IOUint));
    
    // 将 inputParam 的 ioArray 拷贝到 sortedIOArray 中
    for (uint32_t i = 0; i < inputParam->ioVec.len; ++i) {
        sort_IOArray[i] = inputParam->ioVec.ioArray[i];
    }

    // 使用 qsort 对副本进行排序
    qsort(sort_IOArray, inputParam->ioVec.len, sizeof(IOUint), Sort_compare);
    printf("Sort算法得到的io序列");
    // 将排序后的 IOUint 的 id 存入 output->sequence 中
    for (uint32_t i = 0; i < inputParam->ioVec.len; ++i) {
        output->sequence[i] = sort_IOArray[i].id;
        printf("%d ", output->sequence[i]);
    }

    // 释放副本内存
    free(sort_IOArray);
    return 0;
}

// Scan 算法的自定义比较函数：按 startLpos 排序
int Scan_compare(const void *a, const void *b) {
    IOUint *ioA = (IOUint *)a;
    IOUint *ioB = (IOUint *)b;

    if (ioA->startLpos < ioB->startLpos) {
        return -1; // ioA 在 ioB 之前
    } else if (ioA->startLpos > ioB->startLpos) {
        return 1;  // ioA 在 ioB 之后
    } else {
        return 0;  // 两者相等
    }
}

int32_t Scan(const InputParam *inputParam, OutputParam *output){    
    // 创建 ioArray 的副本，用于排序
    IOUint *scan_IOArray = (IOUint *)malloc(inputParam->ioVec.len * sizeof(IOUint));

    // 将 inputParam 的 ioArray 拷贝到 scan_IOArray 中
    for (uint32_t i = 0; i < inputParam->ioVec.len; ++i) {
        scan_IOArray[i] = inputParam->ioVec.ioArray[i];
    }

    // 使用 qsort 对副本进行排序
    qsort(scan_IOArray, inputParam->ioVec.len, sizeof(IOUint), Scan_compare);

    uint32_t seqIndex = 0;
    // 扫描方向从 BOT -> EOT (从小到大)
    for (uint32_t i = 0; i < inputParam->ioVec.len; ++i) {
        if (scan_IOArray[i].wrap % 2 == 0) {  // 假设 wrap 偶数表示从 BOT 向 EOT
            output->sequence[seqIndex++] = scan_IOArray[i].id;
        }
    }

    // 到达 EOT 后，反向扫描（从大到小）
    for (int i = inputParam->ioVec.len - 1; i >= 0; --i) {
        if (scan_IOArray[i].wrap % 2 == 1) {  // 假设 wrap 奇数表示从 EOT 向 BOT
            output->sequence[seqIndex++] = scan_IOArray[i].id;
        }
    }

    printf("Scan算法得到的io序列");
    for (uint32_t i = 0; i < inputParam->ioVec.len; ++i) {
         printf("%d ", output->sequence[i]);
    }

    // 释放副本内存
    free(scan_IOArray);
    return 0;
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
