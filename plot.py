import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 读取Excel文件
df = pd.read_excel('output.xlsx')  # 请将'your_excel_file.xlsx'替换为你的文件名

# 假设第三列为速度列，我们将其命名为 'speed'
speeds = df.iloc[:, 2]  # 获取第三列数据

# 定义速度区间的边界
speed_bins = np.linspace(speeds.min(), speeds.max(), num=30)  # 假设有10个区间，你可以根据需要调整

# 使用numpy的histogram函数计算每个区间的数量
hist, bin_edges = np.histogram(speeds, bins=speed_bins)

# 计算每个区间的中点，用于柱状图的x轴
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# 绘制柱状图
plt.figure(figsize=(10, 6))
plt.bar(bin_centers, hist, width=np.diff(bin_edges), edgecolor='black')

# 设置标题和坐标轴标签
plt.title('Distribution of Speeds')
plt.xlabel('Speed')
plt.ylabel('Number of Occurrences')
#添加图例
plt.legend(['400K'])

# 显示网格
plt.grid(axis='y')

# 显示图形
plt.show()
