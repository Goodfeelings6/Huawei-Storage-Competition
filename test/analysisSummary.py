import os
import re
import pandas as pd

def parse_summary_line(line):
    """ 解析每一行的测试结果 """
    pattern = r"(.+) - (case_\d+\.txt):\s*ioCount:(\d+)\s*algorithmRunningDuration:([0-9.]+)\(ms\)\s*memoryUse:([0-9.]+)\(KB\)\s*addressingDuration:(\d+)\(ms\)\s*readDuration:(\d+)\(ms\)\s*tapeBeltWear:(\d+)\s*tapeMotorWear:(\d+)\s*errorIOCount:(\d+)"
    match = re.match(pattern, line.strip())
    
    if match:
        # 提取并返回每项数据
        exe_name, case_file, io_count, algo_duration, memory_use, addressing_duration, read_duration, tape_belt_wear, tape_motor_wear, error_io_count = match.groups()
        
        return {
            "Executable": exe_name,
            "Case File": case_file,
            "IO Count": int(io_count),
            "Memory Usage (KB)": float(memory_use),
            "Addressing Duration (ms)": int(addressing_duration),
            "Read Duration (ms)": int(read_duration),
            "Tape Belt Wear": int(tape_belt_wear),
            "Tape Motor Wear": int(tape_motor_wear),
            "Error IO Count": int(error_io_count)
        }
    else:
        return None

def read_baseline_summary(baseline_summary_file):
    """ 读取基线程序的 summary.txt，返回每个 case 的性能指标字典 """
    baseline_data = {}
    with open(baseline_summary_file, 'r', encoding='utf-8') as f:
        for line in f:
            result = parse_summary_line(line)
            if result:
                case_file = result["Case File"]
                baseline_data[case_file] = {
                    "Addressing Duration (ms)": result["Addressing Duration (ms)"],
                    "Read Duration (ms)": result["Read Duration (ms)"],
                    "Tape Belt Wear": result["Tape Belt Wear"],
                    "Tape Motor Wear": result["Tape Motor Wear"]
                }
            else:
                print(f"Failed to parse line: {line.strip()}")  # 打印无法解析的行
    print(f"Baseline data loaded: {len(baseline_data)} items")  # 输出基线数据的数量
    return baseline_data



def calculate_percentage_improvement(baseline_value, current_value):
    """ 计算提升百分比 """
    if baseline_value == 0:
        return 0.0  # 防止除以零错误
    return ((baseline_value - current_value) / baseline_value) * 100

def write_to_excel(input_dir, output_excel, baseline_summary_file):
    """ 将 summary.txt 文件的内容和提升百分比写入 Excel 文件 """
    # 读取基线程序的结果
    baseline_data = read_baseline_summary(baseline_summary_file)
    
    result_data = []  # 用来存储所有的测试结果
    
    # 遍历每个可执行程序文件夹中的 A-summary.txt 文件
    for exe_folder in os.listdir(input_dir):
        exe_folder_path = os.path.join(input_dir, exe_folder)
        if os.path.isdir(exe_folder_path):  # 只处理文件夹
            summary_file = os.path.join(exe_folder_path, "A-summary.txt")
            
            if os.path.exists(summary_file):
                with open(summary_file, 'r', encoding='utf-8') as file:
                    for line in file:
                        # 解析每一行数据
                        result = parse_summary_line(line)
                        if result:
                            case_file = result["Case File"]
                            addressing_duration = result["Addressing Duration (ms)"]
                            read_duration = result["Read Duration (ms)"]
                            tape_belt_wear = result["Tape Belt Wear"]
                            tape_motor_wear = result["Tape Motor Wear"]
                            
                            # 查找基线程序的对应值
                            baseline = baseline_data.get(case_file, None)
                            improvement_addressing_duration = improvement_read_duration = improvement_tape_belt_wear = improvement_tape_motor_wear = None
                            
                            if baseline:
                                # 计算各项提升百分比
                                improvement_addressing_duration = calculate_percentage_improvement(baseline["Addressing Duration (ms)"], addressing_duration)
                                improvement_read_duration = calculate_percentage_improvement(baseline["Read Duration (ms)"], read_duration)
                                improvement_tape_belt_wear = calculate_percentage_improvement(baseline["Tape Belt Wear"], tape_belt_wear)
                                improvement_tape_motor_wear = calculate_percentage_improvement(baseline["Tape Motor Wear"], tape_motor_wear)
                            
                            # 更新 result 字典，添加四个新的提升百分比
                            result["Improvement in Addressing Duration (%)"] = improvement_addressing_duration if improvement_addressing_duration is not None else "N/A"
                            result["Improvement in Read Duration (%)"] = improvement_read_duration if improvement_read_duration is not None else "N/A"
                            result["Improvement in Tape Belt Wear (%)"] = improvement_tape_belt_wear if improvement_tape_belt_wear is not None else "N/A"
                            result["Improvement in Tape Motor Wear (%)"] = improvement_tape_motor_wear if improvement_tape_motor_wear is not None else "N/A"
                            
                            # 将结果存储到结果列表中
                            result_data.append(result)
    
    # 使用 pandas 创建 DataFrame 并写入 Excel
    df = pd.DataFrame(result_data)
    
    # 删除不需要的列 'Algorithm Duration (ms)'
    df.drop(columns=["Algorithm Duration (ms)"], inplace=True, errors='ignore')
    
    # 保存为 Excel 文件
    df.to_excel(output_excel, index=False, engine='openpyxl')  # 保存为 Excel 文件
    
    print(f"Data has been written to {output_excel}")

# 调用函数，假设你的输入目录是 args.des，输出文件是 output.xlsx，以及基线程序的 summary.txt
input_directory = '/home/csg/Huawei-Storage-Competition/test/finaltest1210'  # 替换为你实际的目录路径
output_excel = '/home/csg/Huawei-Storage-Competition/test/summaryAnalysisResult.xlsx'  # 输出的 Excel 文件路径
baseline_summary_file = '/home/csg/Huawei-Storage-Competition/test/finaltest1210/project_hw-baseline/A-summary.txt'  # 替换为基线程序的 A-summary.txt 路径

write_to_excel(input_directory, output_excel, baseline_summary_file)

 