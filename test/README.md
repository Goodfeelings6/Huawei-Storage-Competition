# 算例测试

## 目录结构
```
test/
├── README.md 
├── dataset/ 所有测试用例 
├── output/ 默认测试用例输出目录
├── dataset_generator.py 测试用例生成
├── score.py 与基线算法比较计分 脚本
├── test.py 构建或测试指定用例 python脚本
└── testbatch.py 批量测试所有用例 python脚本
```

## 运行测试
### 构建或测试指定用例, test.py 用法：  
```
python test.py s [-src src_path] [-des des_path]
```
+ 参数 s 的可选值为：s0, s1, s2
+ + s0 : 仅构建，即生成可执行文件
+ + s1 : 仅测试，即直接测试指定用例
+ + s2 : 先构建后测试 
+ 参数 -src 为可选参数，用于指定测试用例路径，不指定则默认为 dataset 文件夹下的 case_1.txt 
+ 参数 -des 为可选参数，用于测试输出路径，不指定则默认为 output 文件夹下对应同名文件  

**举例**
```shell
cd test
# 仅构建
python test.py s0
# 仅测试用例 case_1.txt
python test.py s1 -src ./dataset/case_1.txt
# 先构建后测试用例 case_2.txt
python test.py s2 -src ./dataset/case_2.txt
```
---

### 批量测试所有用例, testbatch.py 用法：  
```
python testbatch.py s [-src src_dir] [-des des_dir]
```
+ 参数 s 的可选值为：s0, s1, s2
+ + s0 : 仅构建，即生成可执行文件
+ + s1 : 仅测试，即直接测试所有用例
+ + s2 : 先构建后测试
+ 参数 -src 为可选参数，用于指定测试用例目录，不指定则默认为 dataset 文件夹
+ 参数 -des 为可选参数，用于测试输出目录，不指定则默认为 output 文件夹
+ 注: 测试完成后，还会在输出文件夹输出关键指标汇总文件： A-summary.txt 和关键运行细节汇总文件： A-detail.txt 以及指标和汇总：A-total.txt

**举例**
```shell
cd test
# 仅构建
python testbatch.py s0
# 批量测试所有用例
python testbatch.py s1 
# 先构建后批量测试所有用例
python testbatch.py s2 

# 构建并批量测试所有用例, 输出将位于 output_baseline 文件夹
python testbatch.py s2 -des ./output_baseline
python testbatch.py s2 -des ./output_noshuffle
```
---

### 计算得分, score.py 用法：  
```
python score.py [-c curr_path] [-b base_path]
```
+ 参数 -c 为可选参数，用于指定当前测试算法结果目录，不指定则默认为 output 文件夹
+ 参数 -b 为可选参数，用于指定基线算法结果目录，不指定则默认为 output_baseline 文件夹
+ 注: 输出A-score.txt文件在 -c 参数指定的文件夹

**举例**
```shell
cd test
# 进行计分：即默认对比 output 文件夹 与 output_baseline 文件夹的结果
python score.py 
# 对比 output_other 文件夹 与 output_baseline 文件夹的结果
python score.py -c ./output_other
```
---

### 测试流程
```shell
# 1 修改 algorithm.c 文件，使用 SCAN 基线算法
# 2 每次修改源码后必须重新构建，然后执行批量测试，故使用s2参数，输出文件夹显式指定为 output_baseline
python testbatch.py s2 -des ./output_baseline
# 3 再次修改 algorithm.c 文件，使用 LKH 算法
# 4 同理使用s2参数，输出文件夹采用默认值 output 即可
python testbatch.py s2
# 5 对比计分， 对比的两个结果目录采用默认值即可，可按需求更改
python score.py
```
**注**：如果基线算法不需更改，1，2步 在之后的对比中就可以省略


```shell
# 子集测试
python testbatch.py s2 -src ./dataset_sub -des ./output
python testbatch.py s2 -src ./dataset_sub -des ./output_sub
python score.py -c ./output -b output_sub
```

### 对比批量测试结果, compare.py 用法：  
**更改代码中file1和file2文件名,分别为要对比的两个批量测试结果文件**
**对比结果输出在/home/csg/Huawei-Storage-Competition/test/comparison_output.txt 文件中，更优的结果加双星号表示**