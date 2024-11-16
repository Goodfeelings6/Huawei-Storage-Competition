import matplotlib.pyplot as plt
import os
import re

# 从文件中读取数据
dir = os.path.dirname(__file__)
logFile = os.path.join(dir, 'log.txt')
with open(logFile,'r',encoding='utf-8') as f:
    content = f.read()
# 样例名
name = re.search(r'(case_\d+).txt',content).group(1)
# IO数
ioSize = int(re.search(r'io count = (\d+)',content).group(1))
# 缩放系数值
maxbase = float(re.search(r'maxbase = ([\d.]+)',content).group(1))
# 总时间
totalTime = float(re.search(r'Time.total = ([\d.]+) sec.',content).group(1))
# 预处理时间
preTime = float(re.search(r'Preprocessing time = ([\d.]+) sec.',content).group(1))
# 局部搜索时间
localTime = float(re.search(r'Run 1: Cost = \d+, Time = ([\d.]+) sec.',content).group(1))
# 初始解时间计算
initTime = totalTime-preTime-localTime
# 收敛数据
data = re.findall(r'\* \d+: Cost = (\d+), Time = ([\d.]+) sec.',content)
timeData = []
costData = []
for i,j in data:
    timeData.append(float(j)+initTime+preTime)
    costData.append(int(i))
# 绘图
plt.scatter(timeData, costData, label=f"{name}, IO={ioSize}")
plt.plot(timeData, costData, alpha = 0.5)
plt.text(timeData[0], costData[0], f'begin cost score:\n{round(100-costData[0]*100/maxbase,2)}', fontsize=10, ha='left', va='top')
plt.text(timeData[-1], costData[-1], f'end cost score:\n{round(100-costData[-1]*100/maxbase,2)}', fontsize=10, ha='center', va='bottom')
plt.xlabel("time(s)")
plt.ylabel("cost")
plt.title("Convergence Curve")
plt.legend()
plt.grid(True)  # 显示网格
plt.savefig(os.path.join(dir, "convergence.png"))
plt.show()