#include "algorithm.h"

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

int32_t Sort(const InputParam *input, OutputParam *output){
    // 创建 ioArray 的副本，用于排序
    IOUint *sort_IOArray = (IOUint *)malloc(input->ioVec.len * sizeof(IOUint));
    
    // 将 input 的 ioArray 拷贝到 sortedIOArray 中
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        sort_IOArray[i] = input->ioVec.ioArray[i];
    }

    // 使用 qsort 对副本进行排序
    qsort(sort_IOArray, input->ioVec.len, sizeof(IOUint), Sort_compare);
    // 将排序后的 IOUint 的 id 存入 output->sequence 中
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        output->sequence[i] = sort_IOArray[i].id;
    }
#if DEBUG_LEVEL >= 1
    printf("Sort算法得到的io序列:");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放副本内存
    free(sort_IOArray);
    return 0;
}
