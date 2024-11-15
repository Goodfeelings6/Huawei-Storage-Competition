import argparse
import os
import subprocess
import shutil
import threading
import psutil
import time
import re


# 获取 test 文件夹根路径
testDir = os.path.abspath(os.path.dirname(__file__))
# dataset 文件夹根路径
datasetDir = os.path.join(testDir, "dataset")
# output 文件夹根路径
outputDir = os.path.join(testDir, "output")

# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description="测试单个用例")

# 添加参数
# 必选参数
parser.add_argument("s", 
                    type=str, 
                    choices=["s0", "s1", "s2"],  # 限定参数值范围
                    help="指定仅构建(s0)或直接测试(s1)或先构建后测试(s2)")
# 可选参数，默认值为用例 case_1
parser.add_argument("-src", type=str, help="算例源路径", default=os.path.join(datasetDir, "case_1.txt"))  
# 可选参数，默认值为 outputDir 路径下对应文件
parser.add_argument("-des", type=str, help="算例输出路径", default=None)
# 解析命令行参数
args = parser.parse_args()
# 检查 des 是否使用了默认值
if args.des == None:
    args.__setattr__("des", os.path.join(outputDir, args.src.split('/')[-1]))

memory_usage = []
peak_memory = 0
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
    
def test():
    """ 直接运行测试 """
    if not os.path.exists(os.path.dirname(args.des)): # 若输出文件夹不存在则创建
        os.mkdir(os.path.dirname(args.des))
    exePath = os.path.join(os.path.dirname(testDir), "bin", "project_hw")

    # cmd_str = " ".join([exePath, "-f", args.src, ">"+args.des]) # 使用命令行重定向到对应文件
    cmd_lst = [exePath, "-f", args.src]
    # 启动进程(1若使用shell=True， 故可以直接传递一个命令字符串cmd_str， 但主进程将是sh，由sh派生子进程执行cmd_str命令，
    #         所以统计内存需注意，要统计的是子进程的内存。
    #         2若使用shell=False，只能传递列表命令cmd_lst，此时cmd_lst将作为主进程执行； 其次，无法使用'>'符号来命令行重定向，
    #         只能通过 stdout、 stderr来实现重定向)
    # 为方便统计内存， 此处用 shell=False， stdout重定向
    stdout_file = open(args.des, 'w')
    process = subprocess.Popen(cmd_lst, shell=False, stdout=stdout_file, stderr=subprocess.PIPE, text=True) # 异步的
    
    # 监控内存
    monitor_memory(process, 0.1)
    # 等待进程完成
    process.wait()
    stdout_file.close()
    # 获取进程输出
    stdout, stderr = process.communicate()
    # 检查返回码和输出
    if process.returncode != 0:
        print(f"test {args.src.split('/')[-1]} fail! retval:{process.returncode}, err_detail:{stderr}")
    else:
        print(f"test {args.src.split('/')[-1]} success!")
        # 计算平均内存
        average_memory = sum(memory_usage) / len(memory_usage) if memory_usage else 0
        # print(f"\nAverage Memory Usage: {average_memory:.2f} MB")
        # print(f"Peak Memory Usage: {peak_memory:.2f} MB")
        # 写入对应文件
        # 读取文件内容
        with open(args.des, 'r', encoding='utf-8') as file:
            content = file.read()
        # 替换内存部分内容
        content = re.sub(r"memoryUse:\s+\d+", f"memoryUse:\t\t\t {(average_memory*1024):.2f}", content)
        # 将修改后的内容写回文件
        with open(args.des, 'w', encoding='utf-8') as file:
            file.write(content)

    
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