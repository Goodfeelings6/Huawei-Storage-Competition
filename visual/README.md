## 可视化

### 绘制IO结点的访问顺序图
+ IOGraph.txt 里放入IO算例输入
+ result.txt 里放入对应IO访问序列
+ 执行：
```shell
python visual_sequence.py
```
+ 结果位于同级目录下的 sequence.png 文件


### 绘制算例的收敛曲线
+ log.txt 里放入算例的整个输出日志
+ 执行：
```shell
python visual_convergence.py
```
+ 结果位于同级目录下的 convergence.png 文件