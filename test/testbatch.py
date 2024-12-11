import argparse
import os
import subprocess
import shutil
import re
import threading
import psutil
import time

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
memory_usage = []
peak_memory = 0

# 高性能场景样例hdd
hdd_total_algorithmRunningDuration = 0
hdd_total_addressingDuration = 0
hdd_total_readDuration = 0
hdd_total_readDelay = 0 # 总读时延
hdd_total_tapeBeltWear = 0 # 总带体磨损
hdd_total_tapeMotorWear = 0 # 总电机磨损

# 备分归档场景样例backup
backup_total_algorithmRunningDuration = 0
backup_total_addressingDuration = 0
backup_total_readDuration = 0
backup_total_readDelay = 0 # 总读时延
backup_total_tapeBeltWear = 0
backup_total_tapeMotorWear = 0

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
    # 从ret中顺便提取信息
    global hdd_total_algorithmRunningDuration,hdd_total_addressingDuration,\
           hdd_total_readDuration,hdd_total_tapeBeltWear,hdd_total_tapeMotorWear,\
           backup_total_algorithmRunningDuration,backup_total_addressingDuration,\
           backup_total_readDuration,backup_total_tapeBeltWear,backup_total_tapeMotorWear
    if int(re.search(r'case_(\d+).txt', filePath).group(1))<=60: # 前60个属于hdd场景
        hdd_total_algorithmRunningDuration += float(re.search(r'algorithmRunningDuration:\s*([.\d]+)', ret).group(1))
        hdd_total_addressingDuration += int(re.search(r'addressingDuration:(\d+)\(ms\)', ret).group(1)) 
        hdd_total_readDuration += int(re.search(r'readDuration:(\d+)\(ms\)', ret).group(1)) 
        hdd_total_tapeBeltWear += int(re.search(r'tapeBeltWear:(\d+)', ret).group(1)) 
        hdd_total_tapeMotorWear += int(re.search(r'tapeMotorWear:(\d+)', ret).group(1)) 
    else:
        backup_total_algorithmRunningDuration += float(re.search(r'algorithmRunningDuration:\s*([.\d]+)', ret).group(1))
        backup_total_addressingDuration += int(re.search(r'addressingDuration:(\d+)\(ms\)', ret).group(1)) 
        backup_total_readDuration += int(re.search(r'readDuration:(\d+)\(ms\)', ret).group(1)) 
        backup_total_tapeBeltWear += int(re.search(r'tapeBeltWear:(\d+)', ret).group(1)) 
        backup_total_tapeMotorWear += int(re.search(r'tapeMotorWear:(\d+)', ret).group(1)) 
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
    # match = re.search(r'Subproblem 1 of (\d+)', content)
    # if match:
    #     ret += f"|{match.group(1):^9}"
    # else:
    #     ret += f"|{1:^9}"
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
    # 总用时
    match = re.search(r'algorithmRunningDuration:\s*([\d.]+)\s*\(ms\)', content)
    if match:
        ret += f"|{str(round(float(match.group(1)),2))+'ms':^10}"
    else:
        ret += f"|{'none':^10}ms"
    # 最终Cost
    match = re.search(r'Cost.min\s*=\s*([\d.]+)', content)
    if match:
        ret += f"|{match.group(1):^10}"
    else:
        ret += f"|{'none':^10}"

    ret += "|\n"
    
    global first
    if first: 
        # head = "|   算例   |  IO数量 |子问题数 | Ascent时间|局部搜索最多用时|  最终Cost |\n" \
            #    "|----------|--------|---------|----------|---------------|----------|\n"
        # head = "|   算例   |  IO数量 |  总用时  |  最终Cost |\n" \
        #        "|----------|--------|----------|----------|\n"
        head = "|   算例   |  IO数量 | Ascent时间| 局部搜索 |  总用时  |  最终Cost |\n" \
               "|----------|--------|----------|-----------|----------|----------|\n"
        first = False
        return head+ret
    else:
        return ret

def writeA_total(filePath):
    head = "| 场景 |  总排序时间  |  总寻址时间  | 总读数据时间 |   总读时延   |  总带体磨损  |  总电机磨损  |  加权总消耗    |\n" \
           "|------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|\n"
    with open(filePath, "w", encoding='utf-8') as totalFile:
        totalFile.write(head)
        totalFile.write(f"|{'hdd':^6}")
        totalFile.write(f"|{round(hdd_total_algorithmRunningDuration,2):^14}")
        totalFile.write(f"|{hdd_total_addressingDuration:^14}")
        totalFile.write(f"|{hdd_total_readDuration:^14}")
        hdd_total_readDelay = hdd_total_algorithmRunningDuration+hdd_total_addressingDuration+hdd_total_readDuration
        totalFile.write(f"|{round(hdd_total_readDelay,2):^14}")
        totalFile.write(f"|{hdd_total_tapeBeltWear:^14}")
        totalFile.write(f"|{hdd_total_tapeMotorWear:^14}")
        hdd_weight_sum = 0.5*hdd_total_readDelay+0.3*hdd_total_tapeBeltWear+0.2*hdd_total_tapeMotorWear
        totalFile.write(f"|{round(hdd_weight_sum,2):^14}|\n")

        totalFile.write(f"|{'backup':^6}")
        totalFile.write(f"|{round(backup_total_algorithmRunningDuration,2):^14}")
        totalFile.write(f"|{backup_total_addressingDuration:^14}")
        totalFile.write(f"|{backup_total_readDuration:^14}")
        backup_total_readDelay = backup_total_algorithmRunningDuration+backup_total_addressingDuration+backup_total_readDuration
        totalFile.write(f"|{round(backup_total_readDelay):^14}")
        totalFile.write(f"|{backup_total_tapeBeltWear:^14}")
        totalFile.write(f"|{backup_total_tapeMotorWear:^14}")
        backup_weight_sum = 0.3*backup_total_readDelay+0.5*backup_total_tapeBeltWear+0.2*backup_total_tapeMotorWear
        totalFile.write(f"|{round(backup_weight_sum,2):^14}|\n")

def monitor_memory(process:subprocess.Popen, interval=1):
    """ 监控指定 process 的内存使用 """
    def theadFunc(pid, interval):
        try:
            # 线程通过 pid 获取对应进程信息
            proc = psutil.Process(pid)
            while proc.is_running() and not proc.status() == psutil.STATUS_ZOMBIE:
                # print("主进程：", proc)
                memory_info = proc.memory_info()
                memory_usage_mb = memory_info.rss / 1024 / 1024
                memory_usage.append(memory_usage_mb)
                global peak_memory
                peak_memory = max(peak_memory, memory_usage_mb)
                # print(f"Current Memory: {memory_usage_mb:.2f} MB, Peak Memory: {peak_memory:.2f} MB", end='\r')
                time.sleep(interval)
        except psutil.NoSuchProcess:
            print("Process not found.")

    # 创建一个线程监控 process 的内存使用
    monitor_thread = threading.Thread(target=theadFunc, args=(process.pid, interval))
    monitor_thread.start()

def changeMemoryText(filePath):
    """ 替换输出文件内存部分内容 """
    # 计算平均内存
    average_memory = sum(memory_usage) / len(memory_usage) if memory_usage else 0
    # 读取文件内容
    with open(filePath, 'r', encoding='utf-8') as file:
        content = file.read()
    # 替换内存部分内容
    content = re.sub(r"memoryUse:\s+\d+", f"memoryUse:\t\t\t {(average_memory*1024):.2f}", content)
    # 将修改后的内容写回文件
    with open(filePath, 'w', encoding='utf-8') as file:
        file.write(content)

import os
import subprocess
import re

import os
import subprocess
import re

def test():
    """ 直接运行测试 """
    if not os.path.exists(args.des): # 若输出文件夹不存在则创建
        os.mkdir(args.des)
    
    # 八个可执行文件的路径列表
    exePaths = [
        # "project_hw-nn",
        # "project_hw-gi",
        # "project_hw-ss",
        # "project_hw-mm",
        # "project_hw-lkhnn10",
        # "project_hw-lkhgi10",
        # "project_hw-lkhss10",
        # "project_hw-lkhmm10",
        # "project_hw-lkhnn20",
        # "project_hw-lkhgi20",
        # "project_hw-lkhss20",
        # "project_hw-lkhmm20",        
        # "project_hw-lkhnn30",
        # "project_hw-lkhgi30",
        # "project_hw-lkhss30",
        "project_hw-lkhmm30",
        "project_hw-baseline"
    ]
    
    # 创建每个可执行程序的输出文件夹
    for exe in exePaths:
        exe_dir = os.path.join(args.des, exe)  # 为每个可执行程序创建一个文件夹
        if not os.path.exists(exe_dir):
            os.mkdir(exe_dir)

        # 为每个可执行程序创建独立的 summary, detail 和 total 文件
        summaryFile = open(os.path.join(exe_dir, "A-summary.txt"), "w", encoding='utf-8')
        detailFile = open(os.path.join(exe_dir, "A-detail.txt"), "w", encoding='utf-8')
        totalFile = open(os.path.join(exe_dir, "A-total.txt"), "w", encoding='utf-8')

        # 获取该可执行文件的完整路径
        exePath = os.path.join(os.path.dirname(testDir), "bin", exe)
        print(f"Running tests with executable: {exePath}")

        # 对每个可执行文件进行测试
        for file in sorted(os.listdir(args.src), key=lambda x: int(re.split(r'[_.]+', x)[1])): # 测试所有用例
            # 运行调度算法
            cmd_lst = [exePath, "-f", os.path.join(args.src, file)]
            stdout_file = open(os.path.join(exe_dir, file), 'w')  # 每个文件的输出在对应的文件夹中
            
            # 启动进程
            process = subprocess.Popen(cmd_lst, shell=False, stdout=stdout_file, stderr=subprocess.PIPE, text=True)  # 异步的
            
            # 监控内存
            monitor_memory(process, 0.1)
            # 等待子进程完成
            process.wait()
            stdout_file.close()
            
            # 获取进程输出
            stdout, stderr = process.communicate()

            summaryFile.write(f"{exe} - {file}:\t")
            # 检查返回码和输出
            if process.returncode != 0:
                print(f"test {file} with {exe} fail! retval:{process.returncode}, err_detail:{stderr}")
                summaryFile.write(f"test {file} fail! retval:{process.returncode}, err_detail:{stderr}\n")
                detailFile.write(f"test {file} fail! retval:{process.returncode}, err_detail:{stderr}\n")
            else:
                print(f"test {file} with {exe} success!")
                changeMemoryText(os.path.join(exe_dir, file))
                summaryFile.write(formatMetrics(os.path.join(exe_dir, file)))
                detailFile.write(formatDetails(os.path.join(exe_dir, file)))

            globals()['memory_usage'] = []
            globals()['peak_memory'] = 0

        # 将每个可执行程序的总汇结果写入 A-total.txt
        totalFile.write(f"Results for {exe}:\n")
        totalFile.write(f"Total tests: {len(os.listdir(args.src))}\n")
        # 可以根据需要添加更多汇总信息，或者格式化输出更详细的统计信息
        totalFile.close()

        # 关闭每个程序的 summary 和 detail 文件
        summaryFile.close()
        detailFile.close()

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
