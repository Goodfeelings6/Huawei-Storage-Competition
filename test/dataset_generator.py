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
from functools import cmp_to_key
import os

# 参数, 范围参见 public.h
wrap_size = 280  # 最大wrap规格 280 , [0,278)
tape_size = 730994  # 每个wrap最大lpos规格 730994 , [0,730994)
# 设置随机种子
random.seed(6)
np.random.seed(6)

state_random_probability = 90   # 随机状态触发概率
state_continue_probability = 5 # 连续状态触发概率
state_reverse_probability = 5  # 倒序状态触发概率

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

    # 状态切换
    def state_change():
        state_random = 0   # 随机状态
        state_continue = 0 # 连续状态
        state_reverse = 0  # 倒序状态  
        p = random.randint(0, 100-1)
        if p < state_random_probability: # 触发随机状态
            state_random = 1
        elif p < state_random_probability+state_continue_probability:
            state_continue = random.randint(2, 4) # 触发连续状态,可能2-4个连续
        else:
            state_reverse = random.randint(2, 4) # 触发倒序状态,可能2-4个倒序
        return state_random,state_continue,state_reverse

    state_random,state_continue,state_reverse = 1,0,0 # 初始状态
    # 生成一个 IO 测试用例
    random_io,continue_io,reverse_io = 0,0,0 # 计数
    for i in range(1, io_count + 1):
        if state_random: # 随机状态
            random_io += 1
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
            state_random = state_random - 1
            if state_random == 0:
                state_random,state_continue,state_reverse = state_change()
        elif state_continue: # 连续状态
            continue_io += 1
            first = 1
            while 1:
                io_size = generate_io_size()
                if first: # 定位前一个IO, 只尝试一次，如果失败，将转随机生成
                    wrap = test_cases[-1]['wrap'] 
                    if wrap % 2 == 0: # 正向
                        start_lpos = test_cases[-1]['endLpos']+1
                    else:
                        start_lpos = test_cases[-1]['endLpos']-1
                    first = 0
                else:
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
            state_continue = state_continue - 1
            if state_continue == 0:
                state_random,state_continue,state_reverse = state_change()
        elif state_reverse: # 倒序状态
            reverse_io += 1
            first = 1
            while 1:
                io_size = generate_io_size()
                if first: # 定位前一个IO, 只尝试一次，如果失败，将转随机生成
                    wrap = test_cases[-1]['wrap'] 
                    if wrap % 2 == 0: # 正向
                        end_lpos = test_cases[-1]['startLpos']-1 # 注意倒序状态先确定end_lpos
                    else:
                        end_lpos = test_cases[-1]['startLpos']+1
                    first = 0
                else:
                    wrap = generate_wrap(wrap_distribution)
                    end_lpos = generate_startlpos(startlpos_distribution)
                # 确保 end_lpos 在有效范围内
                if end_lpos < 0:
                    end_lpos = 0
                elif end_lpos > tape_size-1:
                    end_lpos = tape_size-1
                # 计算 start_lpos
                if wrap % 2 == 0: # 正向
                    start_lpos = end_lpos - io_size
                    if start_lpos < 0:  # start_lpos 超出范围
                        start_lpos = 0
                else:
                    start_lpos = end_lpos + io_size
                    if start_lpos > tape_size-1:
                        start_lpos = tape_size-1
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
            state_reverse = state_reverse - 1
            if state_reverse == 0:
                state_random,state_continue,state_reverse = state_change()

    print(f'random_io:{random_io}, continue_io:{continue_io}, reverse_io:{reverse_io}',end='\t')
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

# 排序比较函数
def compare(io1, io2):
    if io1['wrap'] == io2['wrap']:
        if io1['wrap'] % 2 == 0: # 正向
            return io1['startLpos'] - io2['startLpos']
        else:
            return io2['startLpos'] - io1['startLpos']
    else:
        return io1['wrap'] - io2['wrap']


# 测试函数调用
if __name__ == "__main__":
    dataset_dir = './dataset'
    if not os.path.exists(dataset_dir):
        os.mkdir(dataset_dir)

    io_counts = [10,50,100,1000,2000,5000,10000]  # IO 数量
    # io_counts = [1350,1500,1750,3000,4000,6000,7000,8000,9000,9500]  # IO 数量
    wrap_distributions = ['front_20', 'front_30', 'front_50', 'full']  # 在整个磁带区域分布
    startlpos_distributions = ['random', 'gaussian']  # 在横向磁带区域分布

    # 生成测试用例
    index = 1 # 起始命名序号
    # 高性能场景hdd 共60个
    # 随机56个, 带连续的2个, 带连续加倒序的2个
    for io_count in io_counts: # 不同 IO 数量
        state_random_probability = 100   # 随机状态触发概率
        state_continue_probability = 0 # 连续状态触发概率
        state_reverse_probability = 0  # 倒序状态触发概率
        for wrap_distribution in wrap_distributions: # 不同 wrap区域分布
            for startlpos_distribution in startlpos_distributions: # 不同横向磁带区域分布
                # for num in range(5): # 每种生成5个
                for num in range(1): # 每种生成1个
                    test_cases = generate_test_cases(io_count, wrap_distribution, startlpos_distribution)
                    # 写入文件
                    file_name = f'{dataset_dir}/case_{index}.txt'
                    write_test_cases_to_file(test_cases, file_name)
                    print(f'Generate case_{index}.txt success!')
                    index += 1
        # io = 100 和 1000，额外生成带连续的
        if io_count == 100 or io_count == 1000:
            state_random_probability = 80   # 随机状态触发概率
            state_continue_probability = 20 # 连续状态触发概率
            state_reverse_probability = 0  # 倒序状态触发概率
            test_cases = generate_test_cases(io_count, 'full', 'random')
            # 写入文件
            file_name = f'{dataset_dir}/case_{index}.txt'
            write_test_cases_to_file(test_cases, file_name)
            print(f'Generate case_{index}.txt success!')
            index += 1

        # io = 100 和 1000，额外生成带连续加倒序的
        if io_count == 100 or io_count == 1000:
            state_random_probability = 80   # 随机状态触发概率
            state_continue_probability = 10 # 连续状态触发概率
            state_reverse_probability = 10  # 倒序状态触发概率
            test_cases = generate_test_cases(io_count, 'full', 'random')
            # 写入文件
            file_name = f'{dataset_dir}/case_{index}.txt'
            write_test_cases_to_file(test_cases, file_name)
            print(f'Generate case_{index}.txt success!')
            index += 1

    # 备份归档场景backup 共28个
    for io_count in io_counts: # 不同 IO 数量
        state_random_probability = 70   # 随机状态触发概率
        state_continue_probability = 30 # 连续状态触发概率
        state_reverse_probability = 0  # 倒序状态触发概率
        for wrap_distribution in wrap_distributions: # 不同 wrap区域分布
            for num in range(1): # 每种生成1个
                test_cases = generate_test_cases(io_count, wrap_distribution, 'random')
                # 排序
                test_cases = sorted(test_cases, key=cmp_to_key(compare))
                # 随机扰乱 0%~9.9%
                times = int(io_count*random.randint(0, 99)/1000)
                for i in range(times):
                    idx = random.randint(0,io_count-2)
                    # 与后一个io交换
                    test_cases[idx],test_cases[idx+1] = test_cases[idx+1],test_cases[idx]
                # id重排
                for i in range(io_count):
                    test_cases[i]['id'] = i+1
                # 写入文件
                file_name = f'{dataset_dir}/case_{index}.txt'
                write_test_cases_to_file(test_cases, file_name)
                print(f'Generate case_{index}.txt success!')
                index += 1

