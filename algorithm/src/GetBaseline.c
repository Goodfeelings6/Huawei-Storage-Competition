#include "algorithm.h"

// 调用SCAN基线并分析IO序列以获取权重控制系数
void GetBaseline(const InputParam *input, OutputParam *output)
{
    // 获取调度开始时间
    double scheduleStartTime = MyGetTime();
    Scan(input, output);
    /* 基线读时延 */
    base_totaltime += MyGetTime() - scheduleStartTime; // 加上SCAN算法排序时间
    AccessTime accessTime = {0};
    TotalAccessTime(input, output, &accessTime);
    base_totaltime += accessTime.addressDuration; // 加上读寻址时间
    base_totaltime += accessTime.readDuration;    // 加上读数据时间
    /* 基线带体磨损 */
    TapeBeltSegWearInfo segWearInfo = {0};
    base_tapeBeltWear = TotalTapeBeltWearTimes(input, output, &segWearInfo);
    /* 基线电机磨损 */
    base_tapeMotorWear = TotalMotorWearTimes(input, output);

    maxbase = (base_totaltime > base_tapeBeltWear) ? base_totaltime : base_tapeBeltWear;
    maxbase = maxbase * 10;
#if DEBUG_LEVEL >= 1
    printf("base_totaltime = %f\n", base_totaltime);
    printf("base_tapeBeltWear = %f\n", base_tapeBeltWear);
    printf("base_tapeMotorWear = %f\n", base_tapeMotorWear);
    printf("maxbase = %f\n", maxbase);
#endif

    /* 检测算例场景 
       备份归档场景backup：alpha=0.3,beta=0.5,gama=0.2 
       高性能场景hdd：alpha=0.5,beta=0.3,gama=0.2 
    */
    int Nosequential_count=0;//计数非顺序IO节点
    for (uint32_t i = 1; i < input->ioVec.len; ++i)
    {
        if (input->ioVec.ioArray[i].wrap < input->ioVec.ioArray[i - 1].wrap) // 如果后一个wrap小于前一个 肯定非顺序
        {
            Nosequential_count++;
            // printf("Nosequential_id = %d\n", input->ioVec.ioArray[i].id);
        }
        else if (input->ioVec.ioArray[i].wrap == input->ioVec.ioArray[i - 1].wrap)
        {
            if (input->ioVec.ioArray[i].wrap % 2 == 0) // 正向
            {
                if (input->ioVec.ioArray[i].startLpos < input->ioVec.ioArray[i - 1].startLpos)
                    Nosequential_count++;
                    // printf("Nosequential_id = %d\n", input->ioVec.ioArray[i].id);
            }
            else // 反向
            {
                if (input->ioVec.ioArray[i].startLpos > input->ioVec.ioArray[i - 1].startLpos)
                    Nosequential_count++;
                    // printf("Nosequential_id = %d\n", input->ioVec.ioArray[i].id);
            }
        }
        // 是否 hdd 场景
        if (Nosequential_count >= input->ioVec.len * 0.1)
        {
            alpha = 0.5; // 读时延权重
            beta = 0.3;  // 带体磨损权重
            gama = 0.2;  // 电机磨损权重
            break;
        }
    }
    // printf("Nosequential_count = %d\n", Nosequential_count);
#if DEBUG_LEVEL >= 1
    if(alpha==0.3)
        printf("backup\n");
    else
        printf("hdd\n");
#endif
    // 计算最终权重
    final_alpha = alpha / base_totaltime * maxbase;
    final_beta = beta / base_tapeBeltWear * maxbase;
    final_gama = gama / base_tapeMotorWear * maxbase;
#if DEBUG_LEVEL >= 1
    printf("final_alpha = %f\n", final_alpha);
    printf("final_beta = %f\n", final_beta);
    printf("final_gama = %f\n", final_gama);
#endif
}

