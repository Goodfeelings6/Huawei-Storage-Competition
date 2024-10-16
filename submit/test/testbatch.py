import argparse
import os
import subprocess
import shutil
import re


# 获取 test 文件夹根路径
testDir = os.path.abspath(os.path.dirname(__file__))
# dataset 文件夹根路径
datasetDir = os.path.join(testDir, "dataset")
# output 文件夹根路径
outputDir = os.path.join(testDir, "output")

# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description="测试所有用例")
# 必选参数
parser.add_argument("s", 
                    type=str, 
                    choices=["s0", "s1", "s2"],  # 限定参数值范围
                    help="指定仅构建(s0)或直接测试(s1)或先构建后测试(s2)")
# 可选参数，默认值为 dataset 文件夹
parser.add_argument("-src", type=str, help="算例源目录", default=datasetDir)  
# 可选参数，默认值为 output 文件夹
parser.add_argument("-des", type=str, help="算例输出目录", default=outputDir)
# 解析命令行参数
args = parser.parse_args()

first = True


def build():
    """ 仅构建 """
    buildDir = os.path.join(os.path.dirname(testDir), "build")
    if not os.path.exists(buildDir):
        os.mkdir(buildDir)
    else:
        shutil.rmtree(buildDir) # 删除整个build文件夹
        os.mkdir(buildDir)
    # 执行系统命令以构建
    sysCmd = "&&".join([f"cd {buildDir}", "cmake ..", "make"])
    os.system(sysCmd)


def formatMetrics(filePath)->str:
    """ 关键指标输出格式化 """
    ret = ""
    with open(filePath, "r", encoding="utf-8") as f:
        keyFound = 0
        for line in f:
            line = line.strip()
            if line.startswith("Key Metrics"):
                keyFound = 1
            elif keyFound == 0:
                continue
            else: # 读到了关键指标, 如果不想要某个指标, 直接注释对应行就好了
                 if (line.startswith("ioCount")
                    or line.startswith("algorithmRunningDuration")
                    or line.startswith("memoryUse")
                    or line.startswith("addressingDuration")
                    or line.startswith("readDuration")
                    or line.startswith("tapeBeltWear")
                    or line.startswith("tapeMotorWear")
                    or line.startswith("errorIOCount")):
                     line = line.split()
                     ret += "".join(line)
                     ret += "\t"
        ret += "\n"
    return ret

def formatDetails(filePath)->str:
    """ 关键运行细节输出格式化(markdown表格格式) """
    ret = ""
    content = ""
    with open(filePath, "r", encoding="utf-8") as f:
        content = f.read()
    # 算例名
    name = re.split(r'[/.]+',filePath)[-2]
    ret += f"|{name:^10}"
    # IO数量
    match = re.search(r'io count\s*=\s*(\d+)', content)
    if match:
        ret += f"|{match.group(1):^8}"
    else:
        ret += f"|{'none':^8}"
    # 子问题数
    match = re.search(r'Subproblem 1 of (\d+)', content)
    if match:
        ret += f"|{match.group(1):^9}"
    else:
        ret += f"|{1:^9}"
    # Ascent时间
    match = re.search(r'Ascent time\s*=\s*([\d.]+)\s*sec', content)
    if match:
        ret += f"|{match.group(1)+'s':^10}"
    else:
        ret += f"|{'none':^10}s"
    # 局部搜索最多用时
    matchs = re.findall(r'Run 1: Cost = [\d]+, Time = ([\d.]+) sec.', content)
    if matchs:
        maxTime = max(matchs)
        ret += f"|{maxTime+'s':^15}"
    else:
        ret += f"|{'0s':^15}"
    # 最终Cost
    match = re.search(r'addressingDuration:\s*([\d.]+)\s*\(ms\)', content)
    if match:
        ret += f"|{match.group(1):^10}"
    else:
        ret += f"|{'none':^10}"

    ret += "|\n"
    
    global first
    if first: 
        head = "|   算例   |  IO数量 |子问题数 | Ascent时间|局部搜索最多用时|  最终Cost |\n" \
               "|----------|--------|---------|----------|---------------|----------|\n"
        first = False
        return head+ret
    else:
        return ret


def test():
    """ 直接运行测试 """
    if not os.path.exists(args.des): # 若输出文件夹不存在则创建
        os.mkdir(args.des)
    exePath = os.path.join(os.path.dirname(testDir), "build", "project_hw")
    summaryFile = open(os.path.join(args.des, "A-summary.txt"), "w", encoding='utf-8')
    detailFile = open(os.path.join(args.des, "A-detail.txt"), "w", encoding='utf-8')
    for file in sorted(os.listdir(args.src),key=lambda x:int(re.split(r'[_.]+',x)[1])): # 测试所有用例
        # 运行调度算法
        cmd = " ".join([exePath, "-f", os.path.join(args.src, file), ">"+os.path.join(args.des, file)])
        result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
        summaryFile.write(file + ':\t')
        # 检查返回码和输出
        if result.returncode != 0:
            print(f"test {file} fail! retval:{result.returncode}, err_detail:{result.stderr}")
            summaryFile.write(f"test {file} fail! retval:{result.returncode}, err_detail:{result.stderr}\n")
            detailFile.write(f"test {file} fail! retval:{result.returncode}, err_detail:{result.stderr}\n")
        else:
            print(f"test {file} success!")
            summaryFile.write(formatMetrics(os.path.join(args.des, file)))
            detailFile.write(formatDetails(os.path.join(args.des, file)))

    
    summaryFile.close()
    print("Done!")
        
    
def buildAndTest():
    """ 先构建后测试 """
    build()
    test()


if __name__=="__main__":
    if args.s == "s0": 
        build()
    elif args.s == "s1": 
        test()
    elif args.s == "s2": 
        buildAndTest()
