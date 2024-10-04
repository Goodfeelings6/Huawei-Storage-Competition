import argparse
import os
import subprocess
import shutil


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


def test():
    """ 直接运行测试 """
    if not os.path.exists(args.des): # 若输出文件夹不存在则创建
        os.mkdir(args.des)
    exePath = os.path.join(os.path.dirname(testDir), "bin", "project_hw")
    summaryFile = open(os.path.join(args.des, "A-summary.txt"), "w", encoding='utf-8')
    for file in sorted(os.listdir(args.src)): # 测试所有用例
        summaryFile.write(file + ':\t')
        # 运行调度算法
        cmd = " ".join([exePath, "-f", os.path.join(args.src, file), ">"+os.path.join(args.des, file)])
        result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
        # 检查返回码和输出
        if result.returncode != 0:
            print(f"test {file} fail! retval:{result.returncode}, err_detail:{result.stderr}")
            summaryFile.write(f"test {file} fail! retval:{result.returncode}, err_detail:{result.stderr}\n")
        else:
            print(f"test {file} success!")
            summaryFile.write(formatMetrics(os.path.join(args.des, file)))
    
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
