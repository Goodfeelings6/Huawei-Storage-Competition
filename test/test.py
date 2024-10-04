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
    cmd = " ".join([exePath, "-f", args.src, ">"+args.des])
    result = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE, text=True)
    # 检查返回码和输出
    if result.returncode != 0:
        print(f"test {args.src.split('/')[-1]} fail! retval:{result.returncode}, err_detail:{result.stderr}")
    else:
        print(f"test {args.src.split('/')[-1]} success!")
    
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
