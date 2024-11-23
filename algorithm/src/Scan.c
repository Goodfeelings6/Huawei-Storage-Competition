#include "algorithm.h"

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

int32_t Scan(const InputParam *input, OutputParam *output){    
    // 创建 ioArray 的副本，用于排序
    IOUint *scan_IOArray = (IOUint *)malloc(input->ioVec.len * sizeof(IOUint));

    // 将 input 的 ioArray 拷贝到 scan_IOArray 中
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        scan_IOArray[i] = input->ioVec.ioArray[i];
    }

    // 使用 qsort 对副本进行排序
    qsort(scan_IOArray, input->ioVec.len, sizeof(IOUint), Scan_compare);

    uint32_t seqIndex = 0;
    // 扫描方向从 BOT -> EOT (从小到大)
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        if (scan_IOArray[i].wrap % 2 == 0) {  // 假设 wrap 偶数表示从 BOT 向 EOT
            output->sequence[seqIndex++] = scan_IOArray[i].id;
        }
    }

    // 到达 EOT 后，反向扫描（从大到小）
    for (int i = input->ioVec.len - 1; i >= 0; --i) {
        if (scan_IOArray[i].wrap % 2 == 1) {  // 假设 wrap 奇数表示从 EOT 向 BOT
            output->sequence[seqIndex++] = scan_IOArray[i].id;
        }
    }
#if DEBUG_LEVEL >= 1
    printf("Scan算法得到的io序列:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
         printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放副本内存
    free(scan_IOArray);
    return 0;
}
