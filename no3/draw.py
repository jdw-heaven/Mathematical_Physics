import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# 文件路径模板
file_template = '/home/heaven/Desktop/doc/Mathematical_Physics/no3/data/u{}.txt'

# 读取特征值和x坐标
with open('/home/heaven/Desktop/doc/Mathematical_Physics/no3/data/eigenvalue.txt', 'r') as file:
    eigenvalues = list(map(float, file.read().split()))
with open('/home/heaven/Desktop/doc/Mathematical_Physics/no3/data/x.txt', 'r') as file:
    x = list(map(float, file.read().split()))
with open('/home/heaven/Desktop/doc/Mathematical_Physics/no3/data/y.txt', 'r') as file:
    y = list(map(float, file.read().split()))

# 初始化一个字典来存储所有u值
u_values = {}

# 循环读取u1到u16的数据
for i in range(1, 17):  # 从1到16
    f_u = file_template.format(i)
    with open(f_u, 'r') as file:
        u_values[f'u{i}'] = list(map(float, file.read().split()))

# 现在u_values字典包含了所有u值，你可以通过u_values[f'u1']来访问u1的值，以此类推

# 定义子图的数量
num_subplots = 16

# 创建一个新的图形窗口
fig = plt.figure(figsize=(15, 10))  # 可以根据需要调整图形大小

# 随机生成数据点
for i in range(1, num_subplots + 1):
    
    # 添加一个新的子图
    ax = fig.add_subplot(4, 4, i, projection='3d')  # 4行4列的布局，i为当前子图编号
    ax.set_title(f'Eigenvalue {eigenvalues[i-1]}')  # 设置子图标题
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_label(f'u{i}')
    ax.grid(True)
    
    # 画三维散点图
    scatter = ax.scatter(x, y, u_values[f'u{i}'],s=1)

# 调整子图之间的间距
plt.tight_layout()
fig.suptitle(f'Question3,SMR,V = J_0(x_1^(0) rho) ,dots=2601,a=1000.0')
# 显示图形
plt.show()