import re
import pandas as pd

# 从文件中读取内容并解析为字典
def parse_file(filename):
    cases = {}
    with open(filename, 'r', encoding='utf-8') as file:
        for line in file:
            case_name = re.search(r'(case_\d+)', line).group(1)
            addressing_duration = int(re.search(r'addressingDuration:(\d+)', line).group(1))
            cases[case_name] = addressing_duration
    return cases

# 读取两个文件并解析
def compare_addressing_durations(file1, file2, output_file):
    cases1 = parse_file(file1)
    cases2 = parse_file(file2)
    
    # 统计两个文件中较优算例的数量
    file1_better_count = 0
    file2_better_count = 0
    total_duration_file1 = 0
    total_duration_file2 = 0

    # 生成比较结果
    comparison_results = []
    for case in cases1:
        if case in cases2:
            duration1 = cases1[case]
            duration2 = cases2[case]
            total_duration_file1 += duration1
            total_duration_file2 += duration2
            # 比较两个值并标记较优的
            if duration1 < duration2:
                comparison_results.append([f"{case}: **{duration1}**", f"{duration2}"])
                file1_better_count += 1
            elif duration2 < duration1:
                comparison_results.append([f"{case}: {duration1}", f"**{duration2}**"])
                file2_better_count += 1
            else:
                comparison_results.append([f"{case}: {duration1}", f"{duration2}"])

    # 计算两个文件的总差值
    total_difference = total_duration_file1 - total_duration_file2

    # 转换为DataFrame以便写入文件
    df = pd.DataFrame(comparison_results, columns=[f'file1_addressingDuration', f'file2_addressingDuration'])
    
    # 将统计结果和总差值添加到最后两行
    df.loc[len(df.index)] = [f'File 1 better count: {file1_better_count}', f'File 2 better count: {file2_better_count}']
    df.loc[len(df.index)] = [f'Total addressingDuration difference: {total_difference}', '']

    # 将比较结果写入文件
    df.to_csv(output_file, index=False, sep='\t')

# 调用函数进行比较
file1 = '/home/csg/Huawei-Storage-Competition/test/output_SampleSize15P_Trials3/A-summary.txt'
file2 = '/home/csg/Huawei-Storage-Competition/test/output_movetype3-25/A-summary.txt'
output_file = '/home/csg/Huawei-Storage-Competition/test/comparison_output.txt'
compare_addressing_durations(file1, file2, output_file)

print(f"对比结果已生成在 {output_file} 文件中")
