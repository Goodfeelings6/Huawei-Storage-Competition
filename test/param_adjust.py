import os
import re
import subprocess
import shutil

# 样例拷贝，创建子数据集
def copyCase(src, CaseRange:list[tuple], des):
    for i in CaseRange:
        for j in range(i[0],i[1]+1):
            # 样例文件
            file = f"case_{j}.txt"
            # 复制文件
            shutil.copy(os.path.join(src, file), des)  

# 调参信息保存
def saveInfo(scoreFile, infoFile, POPMUSIC_Solutions,POPMUSIC_SampleSize):
    # 读取评分数据
    info = ''
    info += f"|{POPMUSIC_Solutions:^10}"
    info += f"|{POPMUSIC_SampleSize:^8}"
    with open(scoreFile, 'r', encoding='utf-8') as f:
        text = f.read()
        info += "|{:^10}".format(re.search(r'读时延总加分:([\d.]+)',text).group(1))
        info += "|{:^10}".format(re.search(r'带体磨损总加分:([\d.]+)',text).group(1))
        info += "|{:^10}".format(re.search(r'电机磨损总加分:([\d.]+)',text).group(1))
        info += "|{:^10}".format(re.search(r'调度算法总加分:([\d.]+)',text).group(1))
        info += "|{:^8}".format(re.search(r'调度用时总加分:([\d.]+)',text).group(1))
        info += "|{:^8}".format(re.search(r'调度超时总罚分:([\d.]+)',text).group(1))
        info += "|{:^8}".format(re.search(r'空间超限总罚分:([\d.]+)',text).group(1))
        info += "|{:^10}".format(re.search(r'排序错误总数量:([\d.]+)',text).group(1))
        info += "|{:^10}".format(re.search(r'总分:([\d.]+)',text).group(1))
        info += "|\n"

    # 写入文件
    with open(infoFile, 'a', encoding='utf-8') as f:     
        f.write(info)
    
    return float(re.search(r'总分:([\d.]+)',text).group(1)), info

# 调节规模为 ioSize 的算例的某些参数
def adjustParam(ioSize:int, CaseRange:list[tuple], ParamRange:list[tuple]):
    if os.path.exists('./dataset_sub'):
        shutil.rmtree('./dataset_sub')  # 删除整个目录
    os.mkdir('./dataset_sub')
    # 复制对应样例至 ./dataset_sub 文件夹
    copyCase('./dataset', CaseRange, './dataset_sub')

    # 参数循环
    maxScore = 0
    maxInfo = ''
    for POPMUSIC_Solutions in ParamRange[0]:
        for POPMUSIC_SampleSize in ParamRange[1]:
            # 正则匹配替换参数
            with open("../algorithm/algorithm.c",'r+',encoding='utf-8') as f:
                text = f.read()
            newplace = f"//@flag{ioSize}\n\t\tp->POPMUSIC_SampleSize = {POPMUSIC_SampleSize};\n\t\tp->POPMUSIC_Solutions = {POPMUSIC_Solutions};"
            pattern = r"//@flag"+str(ioSize)+r"\s*p->POPMUSIC_SampleSize\s*=\s*\d+;\s*p->POPMUSIC_Solutions\s*=\s*\d+;"
            result = re.sub(pattern, newplace, text)
            # 写回
            with open("../algorithm/algorithm.c",'w',encoding='utf-8') as f:
                f.write(result)

            if os.path.exists('./output_sub'):
                shutil.rmtree('./output_sub')  # 删除整个目录
            # 执行对应批量测试脚本
            cmd = "python testbatch.py s2 -src ./dataset_sub -des ./output_sub"
            process = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
            if process.returncode != 0:
                print(f"testbatch fail! retval:{process.returncode}, err_detail:{process.stderr}")
            
            # 评分
            cmd = "python score.py -c ./output_sub"
            process = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
            if process.returncode != 0:
                print(f"score fail! retval:{process.returncode}, err_detail:{process.stderr}")
            
            # 信息保存
            score, infostr = saveInfo('./output_sub/A-score.txt', './adjustParamInfo.txt',POPMUSIC_Solutions,POPMUSIC_SampleSize)
            if score > maxScore:
                maxScore = score
                maxInfo = infostr

            print(f"io={ioSize},POPMUSIC_Solutions={POPMUSIC_Solutions},POPMUSIC_SampleSize={POPMUSIC_SampleSize} done!")
    # 记录最优
    with open('./adjustParamInfo.txt', 'a', encoding='utf-8') as f:     
        f.write("Optimal:\n")
        f.write(maxInfo)


if __name__=="__main__":
    # 各规模测试样例编号区间(第一个区间为 hdd样例，第二个区间为 backup样例)
    IOCaseRangeDict = { 
        1000:[(27,36),(73,76)],
        2000:[(37,44),(77,80)],
        5000:[(45,52),(81,84)],
        10000:[(53,60),(85,88)]
    }
    # 各规模调参取值(第一个为 POPMUSIC_Solutions 可选取值，第二个区间为 POPMUSIC_SampleSize 可选取值)
    IOParamRangeDict = { 
        1000:[(30,40,50,60),(5,10,15,20)],
        2000:[(30,40,50,60),(5,10,15,20)],
        5000:[(10,15,20,25),(5,10,15,20)],
        10000:[(5,10,15,20),(5,10,15,20)]
    }
    # 要调参的规模
    IOSizes = [1000,2000,5000,10000]
    # IOSizes = [1000]
    # 自动调参
    for ioSize in IOSizes:
        print(f"io={ioSize} start...")
        with open('./adjustParamInfo.txt', 'a', encoding='utf-8') as f:     
            f.write(f"io={ioSize}\n")
            f.write("| solution | sample |读时延加分 | 带体加分 | 电机加分 |调度算法加分|用时加分| 超时罚分| 空间罚分| 排序错误数|  总分   |\n")
            f.write("|----------|--------|----------|----------|----------|----------|--------|--------|--------|----------|----------|\n")
        adjustParam(ioSize, IOCaseRangeDict[ioSize], IOParamRangeDict[ioSize])
    print(f"all done!")
