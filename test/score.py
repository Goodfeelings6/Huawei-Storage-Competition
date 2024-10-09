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
    scoreFile = open(os.path.join(testDir, "score_file.txt"), "w", encoding='utf-8')
    curr_data = []
    base_data = []
    with open(os.path.join(args.c, "A-summary.txt"), "r", encoding="utf-8") as f1:
        curr_data = f1.readlines()
    with open(os.path.join(args.b, "A-summary.txt"), "r", encoding="utf-8") as f2:
        base_data = f2.readlines()


    schedulScore = []           # 调度算法加分
    timeScore = []              # 调度用时加分
    timeoutPenalty = []         # 调度超时罚分
    spaceOverlimitPenalty = []  # 空间超限罚分
    sortErrorPenaly = []        # 排序错误罚分
    for i in range(len(curr_data)): # 比较计算每个算例
        curr_data_i_dict = decode(curr_data[i])
        base_data_i_dict = decode(base_data[i])
        scoreFile.write(curr_data_i_dict["name"] + ':\t')
        # 调度算法加分
        addrT = (base_data_i_dict["addressingDuration"]-curr_data_i_dict["addressingDuration"])*10
        schedulScore.append(addrT)
        # 调度用时加分
        scheduleT = (20000-curr_data_i_dict["algorithmRunningDuration"])*10
        timeScore.append(scheduleT)
        # 

        # 写入文件
        scoreFile.write(f"调度算法加分:{addrT:.2f}\t")
        scoreFile.write(f"调度用时加分:{scheduleT:.2f}\t")
        scoreFile.write('\n')
    
    schedulScore_sum = sum(schedulScore)
    timeScore_sum = sum(timeScore)
    scoreFile.write(f"调度算法总加分:{schedulScore_sum:.2f}\t")
    scoreFile.write(f"调度用时总加分:{timeScore_sum:.2f}\t")

    scoreFile.write(f"\n总分:{schedulScore_sum+timeScore_sum:.2f}\n")
    scoreFile.close()
    print("Done!")
        
if __name__=="__main__":
    score()
