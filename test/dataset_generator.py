""" 
# 测试用例生成脚本
# 规则： 根据数据IO大小、数量、分布特征等因素构造测试用例
1 IO大小——固定大小（ 500lpos、1000lpos、5000lpos ）、随机大小（ 30lpos~1500lpos ）
2 IO数量——10、50、100、1000、2000、5000、10000
3 分布特征
3.1 分布于磁带的区域——分布于磁带前20%、30%、50%、100%
3.2 连续数据段——连续2、3、4个数据段连续
3.3 数据在横向磁带上分布——随机分布、高斯部分
4 边界用例——IO序号连续的数据段、IO序号倒序的数据段
"""

import random
import numpy as np

# 参数, 范围参见 public.h
wrap_size = 280  # 最大wrap规格 280 , [0,278)
tape_size = 730994  # 每个wrap最大lpos规格 730994 , [0,730994)
# 设置随机种子
random.seed(6)
np.random.seed(6)

# 测试用例生成函数
def generate_test_cases(io_count, wrap_distribution, startlpos_distribution):
    test_cases = []
    wrap_flag_dic = {} # 存每个wrap已选区间列表
    for i in range(wrap_size):
        wrap_flag_dic[i] = []
    
    # IO 大小生成
    def generate_io_size():
        if random.randint(0, 100) < 15:
            # 15% 几率固定大小
            return random.choice([500, 1000, 5000])
        else:
            # 否则随机大小
            return random.randint(30, 1500)

    # wrap 生成，给定分布于磁带的哪些区域
    def generate_wrap(distribution): 
        if distribution == 'front_20':
            return random.randint(0, int((wrap_size-1) * 0.2))
        elif distribution == 'front_30':
            return random.randint(0, int((wrap_size-1) * 0.3))
        elif distribution == 'front_50':
            return random.randint(0, int((wrap_size-1) * 0.5))
        elif distribution == 'full':
            return random.randint(0, (wrap_size-1))
    
    # startlpos 生成，给定横向磁带上的分布
    def generate_startlpos(distribution):
        if distribution == 'random':
            return random.randint(0, tape_size-1)
        elif distribution == 'gaussian':
            # 以横向磁带中间为中心的高斯分布
            return int(np.random.normal(loc=(tape_size-1)//2, scale=(tape_size-1)//6))
        
    # 生成一个 IO 测试用例
    for i in range(1, io_count + 1):
        while 1:
            io_size = generate_io_size()
            wrap = generate_wrap(wrap_distribution)
            start_lpos = generate_startlpos(startlpos_distribution)
            # 确保 start_lpos 在有效范围内
            if start_lpos < 0:
                start_lpos = 0
            elif start_lpos > tape_size-1:
                start_lpos = tape_size-1
            # 计算 end_lpos
            if wrap % 2 == 0: # 正向
                end_lpos = start_lpos + io_size
                if end_lpos > tape_size-1:  # end_lpos 超出范围
                    end_lpos = tape_size-1
            else:
                end_lpos = start_lpos - io_size
                if end_lpos < 0:
                    end_lpos = 0
            # 如果起始与结束重合，跳过
            if start_lpos == end_lpos:
                continue
            
             # 检查新生成的子区间是否与之前的子区间重叠
            seg = (start_lpos, end_lpos)
            if not any(start_lpos < existing_end and end_lpos > existing_start 
                       for existing_start, existing_end in wrap_flag_dic[wrap]):
                # 区间未重叠， 可行
                test_cases.append({
                    'id': i,
                    'wrap': wrap,
                    'startLpos': start_lpos,
                    'endLpos': end_lpos
                })
                wrap_flag_dic[wrap].append(seg) # 记录这个区间
                break

    return test_cases

# 随机生成磁头信息
def generate_head_info():
    # 随机生成 head 所处 wrap 值 
    head = random.randint(0, wrap_size-1)
    # 随机生成 lpos 值
    lpos = random.randint(0, tape_size-1) 
    # 随机生成 status 值（0 静止，1 读写）
    status = random.randint(0, 1)
    return head, lpos, status

# 将测试用例写入文件
def write_test_cases_to_file(test_cases, file_name):
    with open(file_name, 'w') as f:
        # 随机生成并写入磁头信息
        head, lpos, status = generate_head_info()
        f.write(f'["head":"wrap","lpos","status"]\n')
        f.write(f'[{head},{lpos},{status}]\n')
        
        # 写入 IO 数量
        f.write('["io count"]\n')
        f.write(f'[{len(test_cases)}]\n')
        
        # 写入 IO 信息
        f.write('["io":"id","wrap","startLpos","endLpos"]\n')
        for case in test_cases:
            f.write(f'[{case["id"]},{case["wrap"]},{case["startLpos"]},{case["endLpos"]}]\n')

# 测试函数调用
if __name__ == "__main__":
    io_counts = [10,50,100,1000,2000,5000,10000]  # IO 数量
    # io_counts = [1350,1500,1750,3000,4000,6000,7000,8000,9000,9500]  # IO 数量
    wrap_distributions = ['front_20', 'front_30', 'front_50', 'full']  # 在整个磁带区域分布
    startlpos_distributions = ['random', 'gaussian']  # 在横向磁带区域分布

    # 生成测试用例
    index = 1 # 起始命名序号
    # index = 286 # 起始命名序号
    for io_count in io_counts: # 不同 IO 数量
        for wrap_distribution in wrap_distributions: # 不同 wrap区域分布
            for startlpos_distribution in startlpos_distributions: # 不同横向磁带区域分布
                # for num in range(5): # 每种生成5个
                for num in range(1): # 每种生成1个
                    test_cases = generate_test_cases(io_count, wrap_distribution, startlpos_distribution)
                    # 写入文件
                    file_name = f'./dataset/case_{index}.txt'
                    write_test_cases_to_file(test_cases, file_name)
                    print(f'Generate case_{index}.txt success!')
                    index += 1

