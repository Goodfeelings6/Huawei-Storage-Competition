import argparse
import os
import re


# 获取 test 文件夹根路径
testDir = os.path.abspath(os.path.dirname(__file__))
# output 文件夹根路径
outputDir = os.path.join(testDir, "output")
# output_baseline 文件夹根路径
outputBaseDir = os.path.join(testDir, "output_baseline")

# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description="计算得分")

# 可选参数，默认值为 dataset 文件夹
parser.add_argument("-c", type=str, help="当前测试算法的结果目录", default=outputDir)  
# 可选参数，默认值为 output 文件夹
parser.add_argument("-b", type=str, help="基线算法的结果目录", default=outputBaseDir)
# 解析命令行参数
args = parser.parse_args()

def decode(s:str)->dict:
    """ 关键指标解码为字典格式 """
    ret = {}
    split_s = [item for item in re.split(r'[\(\)\s:]+', s) if item]
    # print(split_s)
    ret["name"] = split_s[0]
    i = 1
    while i < len(split_s):
        if split_s[i] == 'ms' or split_s[i] == 'KB':
            i += 1
            continue
        else:
            ret[split_s[i]] = float(split_s[i+1])
            i += 2
    return ret


def score():
    """ 计算得分 """
    scoreFile = open(os.path.join(args.c, "A-score.txt"), "w", encoding='utf-8')
    curr_data = []
    base_data = []
    with open(os.path.join(args.c, "A-summary.txt"), "r", encoding="utf-8") as f1:
        curr_data = f1.readlines()
    with open(os.path.join(args.b, "A-summary.txt"), "r", encoding="utf-8") as f2:
        base_data = f2.readlines()

    readTimeScoreList = []      # 读时延加分
    tapeBeltWearScoreList = []  # 带体磨损加分
    tapeMotorWearScoreList = [] # 电机磨损加分
    schedulScore = []           # 调度算法加分
    timeScore = []              # 调度用时加分
    timeoutPenalty = []         # 调度超时罚分
    spaceOverlimitPenalty = []  # 空间超限罚分
    sortErrorCount = []         # 排序错误数量
    instanceScore = []          # 算例总得分
    for i in range(len(curr_data)): # 比较计算每个算例
        curr_data_i_dict = decode(curr_data[i])
        instanceID = int(re.search(r'case_(\d+).txt', curr_data_i_dict["name"]).group(1))
        base_data_i_dict = decode(base_data[instanceID-1])
        # 权重
        if instanceID <= 60:
            alpha = 0.5
            beta = 0.3
            gama = 0.2
        else:
            alpha = 0.3
            beta = 0.5
            gama = 0.2
        scoreFile.write(curr_data_i_dict["name"] + ':\t')
        scoreFile.write("IO:"+str(curr_data_i_dict["ioCount"])+"\t")
        # 调度算法加分
        base_read_time = base_data_i_dict["algorithmRunningDuration"]+base_data_i_dict["addressingDuration"]+base_data_i_dict["readDuration"]
        curr_read_time = curr_data_i_dict["algorithmRunningDuration"]+curr_data_i_dict["addressingDuration"]+curr_data_i_dict["readDuration"]
        ## 三部分得分
        readTimeScore = alpha*(base_read_time-curr_read_time)/base_read_time*100
        tapeBeltWearScore = beta*(base_data_i_dict["tapeBeltWear"]-curr_data_i_dict["tapeBeltWear"])/base_data_i_dict["tapeBeltWear"]*100
        tapeMotorWearScore = gama*(base_data_i_dict["tapeMotorWear"]-curr_data_i_dict["tapeMotorWear"])/base_data_i_dict["tapeMotorWear"]*100
        readTimeScoreList.append(readTimeScore)
        tapeBeltWearScoreList.append(tapeBeltWearScore)
        tapeMotorWearScoreList.append(tapeMotorWearScore)
        ## 相加
        addrT = readTimeScore+tapeBeltWearScore+tapeMotorWearScore
        if addrT < 0 :
            addrT = 0
        schedulScore.append(addrT)
        # 调度用时加分 和 调度超时罚分
        if curr_data_i_dict["algorithmRunningDuration"]<=20000:
            scheduleT = (20000-curr_data_i_dict["algorithmRunningDuration"])/20000*10*(curr_data_i_dict["ioCount"]/10000)
            penaltyT = 0
        else:
            scheduleT = 0
            penaltyT = (curr_data_i_dict["algorithmRunningDuration"]-20000)/20000*10*(2-curr_data_i_dict["ioCount"]/10000)
        timeScore.append(scheduleT)
        timeoutPenalty.append(penaltyT)
        # 空间超限罚分
        if curr_data_i_dict["memoryUse"] <= 10240:
            spaceOverlimitP = 0
        else:
            spaceOverlimitP = (curr_data_i_dict["memoryUse"]-10240)/10240*10*(2-curr_data_i_dict["ioCount"]/10000)
        spaceOverlimitPenalty.append(spaceOverlimitP)
        # 排序错误数量
        errorC = curr_data_i_dict["errorIOCount"]
        sortErrorCount.append(errorC)

        # 算例得分
        if errorC > 0:
            instanceS = 0
        else:
            instanceS = addrT+scheduleT-penaltyT-spaceOverlimitP
        instanceScore.append(instanceS)

        # 写入文件
        scoreFile.write(f"读时延加分:{readTimeScore:.2f}\t")
        scoreFile.write(f"带体磨损加分:{tapeBeltWearScore:.2f}\t")
        scoreFile.write(f"电机磨损加分:{tapeMotorWearScore:.2f}\t")
        scoreFile.write(f"调度算法加分:{addrT:.2f}\t")
        scoreFile.write(f"调度用时加分:{scheduleT:.2f}\t")
        scoreFile.write(f"调度超时罚分:{penaltyT:.2f}\t")
        scoreFile.write(f"空间超限罚分:{spaceOverlimitP:.2f}\t")
        scoreFile.write(f"排序错误数量:{errorC:.2f}\t")
        scoreFile.write(f"算例总得分:{instanceS:.2f}\t")
        scoreFile.write('\n')
    
    scoreFile.write(f"读时延总加分:{sum(readTimeScoreList):.2f}\t")
    scoreFile.write(f"带体磨损总加分:{sum(tapeBeltWearScoreList):.2f}\t")
    scoreFile.write(f"电机磨损总加分:{sum(tapeMotorWearScoreList):.2f}\t")
    scoreFile.write(f"调度算法总加分:{sum(schedulScore):.2f}\t")
    scoreFile.write(f"调度用时总加分:{sum(timeScore):.2f}\t")
    scoreFile.write(f"调度超时总罚分:{sum(timeoutPenalty):.2f}\t")
    scoreFile.write(f"空间超限总罚分:{sum(spaceOverlimitPenalty):.2f}\t")
    scoreFile.write(f"排序错误总数量:{sum(sortErrorCount):.2f}\t")
    scoreFile.write(f"\n总分:{sum(instanceScore):.2f}\n")

    # 测相比基线的提升幅度
    # with open(os.path.join(args.c, "A-total.txt"), "r", encoding="utf-8") as f1:
    #     curr_total_cost = float((f1.readlines()[-1]).strip().split('|')[-2].strip())
    # with open(os.path.join(args.b, "A-total.txt"), "r", encoding="utf-8") as f2:
    #     base_total_cost = float((f2.readlines()[-1]).strip().split('|')[-2].strip())
    # scoreFile.write(f"优化幅度:{(base_total_cost-curr_total_cost)/base_total_cost*100:.3f}%\n")
    scoreFile.close()
    print("Score Done!")
        
if __name__=="__main__":
    score()
