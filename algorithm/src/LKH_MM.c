#include "algorithm.h"

// LKH 算法 + 最大匹配(MM)生成初始解
int32_t LKH_MM(const InputParam *input, OutputParam *output)
{   
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();

    int32_t ret = 0;

    uint32_t len = output->len + 2; //增加了一个虚拟节点和磁头起始位置节点

    // 最大匹配获取初始解
    MaxMatching_Greedy(input, output);
    // 路径拷贝至初始解数组
    int *intialTour = (int *)malloc(len * sizeof(int));
    intialTour[0] = len - 1, intialTour[1] = len;
    for (int i = 2; i < len; i++)
        intialTour[i] = output->sequence[i - 2];

    /* 设置固定边, 结点从1开始编号 */
    int fixLen = 1;
    int *fixEdge = (int *)malloc(2 * fixLen * sizeof(int));
    // 虚拟结点到磁头结点的边固定
    fixEdge[0] = len - 1;
    fixEdge[1] = len;

    /* 调用 LKH 求解 */
    /* 确定LKH输入结构体 */
    LKHInput *lkhInput = (LKHInput *)malloc(sizeof(LKHInput));
    lkhInput->adjMat = 0;
    lkhInput->matDimension = len;
    lkhInput->intialTour = intialTour;
    lkhInput->fixEdge = fixEdge;
    lkhInput->fixEdgeLen = fixLen;
    lkhInput->scheduleStartTime = scheduleStartTime;
    lkhInput->lkhParameters = (LKHParameters *)malloc(sizeof(LKHParameters));
    loadDefaultParam(lkhInput->lkhParameters);
    loadUserChangedParam(lkhInput->matDimension, lkhInput->lkhParameters, scheduleStartTime);
    lkhInput->lkhParameters->OriginInput = input;
    /* 确定LKH输出结构体 */
    LKHOutput *lkhOutput = (LKHOutput *)malloc(sizeof(LKHOutput));
    lkhOutput->tourCost = 0;
    lkhOutput->tourResult = (int *)malloc(len * sizeof(int));

    ret = solveTSP(lkhInput, lkhOutput);
   
    // 处理解
    int i, j;
    for (i = 0; i < len; ++i){
        if(lkhOutput->tourResult[i] == len) // i指向磁头结点
            break;
    }
    for (j = 0; j < len - 2; ++j){ // 排序结果赋值到输出
        if(i + 1 > len - 1)
            i = -1;
        output->sequence[j] = lkhOutput->tourResult[++i];
    }

    /* 打印求解结果：输出output->sequence的内容 */
#if DEBUG_LEVEL >= 1
    printf("LKH result\n");
    printf("io len:%d\n",len-2);
    printf("Output sequence:\n");
    for (uint32_t i = 0; i < len - 2; ++i) {
        printf("%d ", output->sequence[i]);
    }
    printf("\n");
#endif

    free(intialTour);
    free(fixEdge);
    free(lkhInput->lkhParameters);
    free(lkhInput);
    free(lkhOutput->tourResult);
    free(lkhOutput);

    return ret;
}
