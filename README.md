## 开发环境要求：
Ubuntu  ≥ 	18.04
cmake	≥  	3.27.5
make	≥ 	4.1
GCC		≥ 	7.5

## 编译运行

```shell
# 命令
cd project_hw
mkdir build
cd build
cmake ..
make
../bin/project_hw -f ../test/dataset/case_1.txt
```

```shell
cd build
rm -rf *
cmake ..
make
../bin/project_hw -f ../test/dataset/case_1.txt
```

## 测试
### 先入先出算法
| 算例      | 寻址时间(ms) |
|----------|--------------|
| case1    |  1362738     |
| case2    |  850195      |
| case3    |  1003042     |
| case4    |  3859545     |
| case5    |  7664161     |

### 无虚拟结点ATSP, LKH2, 默认值
| 算例      | 寻址时间(ms) |
|----------|--------------|
| case1    |  623643      |
| case2    |  429979      |
| case3    |  545996      |
| case4    |  676173      |
| case5    |  1082234     |

### 虚拟结点ATSP, LKH2, 默认值
| 算例      | 寻址时间(ms) |
|----------|--------------|
| case1    |  354295      |
| case2    |  293696      |
| case3    |  495332      |
| case4    |  685706      |
| case5    |  996508      |