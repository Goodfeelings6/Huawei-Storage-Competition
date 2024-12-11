#include "algorithm.h"

// 按 startLpos 排序
int SS_compare(const void *a, const void *b) {
    IOUintSS *ioA = (IOUintSS *)a;
    IOUintSS *ioB = (IOUintSS *)b;

    if (ioA->startLpos < ioB->startLpos) {
        return -1; // ioA 在 ioB 之前
    } else if (ioA->startLpos >= ioB->startLpos) {
        return 1;  // ioA 在 ioB 之后
    } else {
        return 0;  // 两者相等
    }
}

// 多轮线性SCAN构造-线性搜索
int32_t MultiSortScan(const InputParam *input, OutputParam *output){    
    // 统计正反wrap上IO数量
    int oddWrapIOSize=0; // 反向wrap(奇数)上的IO数量
    int evenWrapIOSize=0; // 正向wrap(偶数)上的IO数量
    for (uint32_t i = 0; i < input->ioVec.len; ++i) 
    {
        if(input->ioVec.ioArray[i].wrap%2==1)oddWrapIOSize++;
        else evenWrapIOSize++;
    }
    // 创建 ioArray 的副本，用于排序
    IOUintSS *odd_IOArray = (IOUintSS *)malloc(oddWrapIOSize* sizeof(IOUintSS));
    IOUintSS *even_IOArray  = (IOUintSS *)malloc(evenWrapIOSize * sizeof(IOUintSS));

    uint32_t oddIndex = 0;
    uint32_t evenIndex = 0;
    
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
        if(input->ioVec.ioArray[i].wrap%2==1)
        {
            odd_IOArray[oddIndex].endLpos = input->ioVec.ioArray[i].endLpos+1;
            odd_IOArray[oddIndex].startLpos = input->ioVec.ioArray[i].startLpos+1;
            odd_IOArray[oddIndex].wrap=input->ioVec.ioArray[i].wrap;
            odd_IOArray[oddIndex].id=input->ioVec.ioArray[i].id;
            odd_IOArray[oddIndex].visit=0;
            oddIndex++;
        }
        else
        {
            even_IOArray[evenIndex].endLpos = input->ioVec.ioArray[i].endLpos+1;
            even_IOArray[evenIndex].startLpos = input->ioVec.ioArray[i].startLpos+1;
            even_IOArray[evenIndex].wrap=input->ioVec.ioArray[i].wrap;
            even_IOArray[evenIndex].id=input->ioVec.ioArray[i].id;
            even_IOArray[evenIndex].visit=0;
            evenIndex++;
        }
    }
    
    // 使用 qsort 对副本进行排序
    qsort(odd_IOArray, oddIndex, sizeof(IOUintSS), SS_compare);
    qsort(even_IOArray, evenIndex, sizeof(IOUintSS), SS_compare);
    
    // 多轮SortScan扫描构造
    uint32_t seqIndex = 0;
    int current_wrap = input->headInfo.wrap; // 当前所处wrap
    int current_lpos = input->headInfo.lpos+1; // 当前所处lpos

    while(seqIndex<input->ioVec.len)
    {  
        // 正向扫描
        for (uint32_t i = 0; i < evenWrapIOSize; ++i)
        { 
            if((even_IOArray[i].visit==0)&&(even_IOArray[i].startLpos>current_lpos))
            {   
                // 如果存在多个startLpos相同的IO，选择离当前所处wrap最近的
                int index=i+1;
                while(index<=evenWrapIOSize-1&&even_IOArray[index].startLpos==even_IOArray[i].startLpos)
                {     
                    if((even_IOArray[index].visit==0)&&abs((int)even_IOArray[index].wrap-current_wrap)<abs((int)even_IOArray[i].wrap-current_wrap))
                    {
                        i=index;
                    }
                    index++;
                }
                even_IOArray[i].visit=1;
                current_wrap=even_IOArray[i].wrap;
                current_lpos=even_IOArray[i].endLpos;
                output->sequence[seqIndex++]=even_IOArray[i].id;
            }
        }
        // 转反向扫描
        current_lpos = odd_IOArray[oddWrapIOSize-1].startLpos+1;
        current_wrap = odd_IOArray[oddWrapIOSize-1].wrap;
        for (int i = oddWrapIOSize-1; i >= 0; i--)
        {   
            if(odd_IOArray[i].visit==0&&odd_IOArray[i].startLpos<current_lpos)
            {   
                // 如果存在多个startLpos相同的IO，选择离当前所处wrap最近的
                int index=i-1;
                while((index>=0)&&odd_IOArray[index].startLpos==odd_IOArray[i].startLpos)
                {     
                    if((odd_IOArray[index].visit==0)&&abs((int)odd_IOArray[index].wrap-current_wrap)<abs((int)odd_IOArray[i].wrap-current_wrap))
                    {
                        i=index;
                    }
                    index--;
                }
                odd_IOArray[i].visit=1;
                current_wrap=odd_IOArray[i].wrap;
                current_lpos=odd_IOArray[i].endLpos;
                output->sequence[seqIndex++]=odd_IOArray[i].id;
            }
        }
        // 转正向扫描
        current_lpos = even_IOArray[0].startLpos-1;
        current_wrap = even_IOArray[0].wrap;
    }
    
#if DEBUG_LEVEL >= 1
    printf("多轮线性SCAN构造-线性搜索:\n");
    for (uint32_t i = 0; i < input->ioVec.len; ++i) {
         printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    // 释放副本内存
    free(odd_IOArray);
    free(even_IOArray);
    return 0;
}
