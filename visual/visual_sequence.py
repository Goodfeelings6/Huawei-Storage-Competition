import matplotlib.pyplot as plt
import os

# 读取数据
def read_io_data(filename):
    io_requests = []
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i in range(5, len(lines)):  # 从第5行开始是IO请求数据
            if lines[i].startswith('["io"'):
                continue  # 跳过表头
            parts = lines[i].strip().strip('[').strip(']').split(',')
            io_requests.append({
                "id": int(parts[0]),
                "wrap": int(parts[1]),
                "start": int(parts[2]),
                "end": int(parts[3])
            })
    return io_requests

# 读取 IO 请求序列（格式为: 8 2 5 7 1 4 9 3 10 6）
def read_io_sequence(filename):
    with open(filename, 'r') as file:
        sequence = [int(x) for x in file.read().split()]
    return sequence

# 从文件中读取IO请求数据
dir = os.path.dirname(__file__)
io_requests_with_wrap = read_io_data(os.path.join(dir, 'IOGraph.txt'))  
io_sequence = read_io_sequence(os.path.join(dir,'result.txt'))  # 读取指定的IO序列

# 提取wrap请求数据
wrap_requests = [(req["id"], req["wrap"], req["start"], req["end"]) for req in io_requests_with_wrap]

# 创建图表并根据wrap奇偶性显示不同颜色的IO请求
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

# 提取所有wrap行
wraps = sorted(set([req['wrap'] for req in io_requests_with_wrap]))

# 绘制每个wrap行的横线
for wrap in wraps:
    ax.axhline(y=wrap, color='gray', linestyle='--', alpha=0.5)

# 绘制每个请求的块范围，奇数wrap使用蓝色，偶数wrap使用红色
for req in wrap_requests:
    id_, wrap, start, end = req
    color = 'blue' if wrap % 2 == 1 else 'red'
    ax.plot([start, end], [wrap, wrap], marker='o', color=color, label=f"IO {id_} (Wrap {wrap})")

# 获取每个IO请求的中心点
def get_io_center(request):
    return (request["start"] + request["end"]) / 2, request["wrap"]

# 按序连线
sequence_points = []
for io_id in io_sequence:
    # 查找该IO请求
    request = next(req for req in io_requests_with_wrap if req["id"] == io_id)
    # 获取中心点
    center_x, center_y = get_io_center(request)
    sequence_points.append((center_x, center_y))

# 将中心点连线
for i in range(1, len(sequence_points)):
    x_values = [sequence_points[i-1][0], sequence_points[i][0]]
    y_values = [sequence_points[i-1][1], sequence_points[i][1]]
    ax.plot(x_values, y_values, color='green', linestyle='-', marker='x')

# 设置图表标题和标签

ax.set_xlabel('lpos')
ax.set_ylabel('Wrap Line')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig(os.path.join(dir, "sequence.png"))
plt.show()
