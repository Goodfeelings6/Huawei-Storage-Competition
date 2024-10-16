import os

def merge_files(input_dir, output_file):
    # 获取输入目录下的所有文件
    files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
    files.sort(key=lambda x:x.lower())
    with open(output_file, 'w') as outfile:
        for filename in files:
            file_path = os.path.join(input_dir, filename)
            # 添加文件名标题
            outfile.write(f"//################################### {filename} begin ###################################\n")
            # 读取并写入文件内容
            with open(file_path, 'r') as infile:
                outfile.write(infile.read())
            # 添加文件结束标识（可选）
            outfile.write(f"//################################### {filename} end ###################################\n\n")

if __name__ == "__main__":
    input_directory = "/home/csg/Huawei-Storage-Competition/LKH2_MIN/src"  # 修改为你的输入文件夹路径
    output_filename = "merged_output.txt"  # 输出合并后的文件名
    merge_files(input_directory, output_filename)
    print(f"All files merged into {output_filename}")
