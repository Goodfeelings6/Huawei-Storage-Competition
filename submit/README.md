# “Massive Storage”第二届大学生信息存储技术竞赛·挑战赛

## 开发环境要求：
Ubuntu  ≥ 	18.04
cmake	≥  	3.27.5
make	≥ 	4.1
GCC		≥ 	7.5

## 编译运行

```shell
# 命令
mkdir build
cd build
cmake ..
make
./project_hw -f ../dataset/case_1.txt
```


## 批量测试 

```shell
cd test
# 构建并批量测试 ../dataset 内所有用例, 输出将位于 output 文件夹
python testbatch.py s2 -src ../dataset -des ./output
```