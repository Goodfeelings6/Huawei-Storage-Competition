import psutil
import time
import subprocess

def monitor_memory(command):
    memory_usage = []
    peak_memory = 0

    try:
        # 启动进程
        process = subprocess.Popen(command, shell=True)
        pid = process.pid
        proc = psutil.Process(pid)

        print(f"Monitoring memory usage for process PID: {pid}")
        
        # 持续监控内存
        while process.poll() is None:  # 检查进程是否结束
            try:
                mem_info = proc.memory_info()
                current_memory = mem_info.rss / (1024 ** 2)  # 转换为 MB
                memory_usage.append(current_memory)
                peak_memory = max(peak_memory, current_memory)

                print(f"Current Memory: {current_memory:.2f} MB, Peak Memory: {peak_memory:.2f} MB", end='\r')
                time.sleep(1)  # 每秒更新一次
            except psutil.NoSuchProcess:
                #print("\nProcess has terminated.")
                break

        # 计算平均内存
        average_memory = sum(memory_usage) / len(memory_usage) if memory_usage else 0
        print(f"\nAverage Memory Usage: {average_memory:.2f} MB")
        print(f"Peak Memory Usage: {peak_memory:.2f} MB")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # 指定执行命令
    command = "../bin/project_hw -f ../test/dataset/case_50.txt"
    monitor_memory(command)
