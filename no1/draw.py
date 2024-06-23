import matplotlib.pyplot as plt
import numpy as np

# 文件路径
f_eigenvalue = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/eigenvalue.txt'
f_x = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/x.txt'
f_u1 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u1.txt'
f_u2 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u2.txt'
f_u3 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u3.txt'
f_u4 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u4.txt'
f_u5 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u5.txt'
f_u6 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u6.txt'
f_u7 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u7.txt'
f_u8 = '/home/heaven/Desktop/doc/Mathematical_Physics/no1/data/u8.txt'

# 读取数据并转换为浮点数列表
with open(f_eigenvalue, 'r') as file:
    eigenvalue = list(map(float, file.read().split()))  # 假设只有一个特征值
with open(f_x, 'r') as file:
    x = list(map(float, file.read().split()))  # 读取x坐标并转换为浮点数
with open(f_u1, 'r') as file:
    u1 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u2, 'r') as file:
    u2 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u3, 'r') as file:
    u3 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u4, 'r') as file:
    u4 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u5, 'r') as file:
    u5 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u6, 'r') as file:
    u6 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u7, 'r') as file:
    u7 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数
with open(f_u8, 'r') as file:
    u8 = list(map(float, file.read().split()))  # 读取u值并转换为浮点数

# 创建一个2行4列的子图布局
fig, axs = plt.subplots(2, 4)

axs[0, 0].scatter(x, u1,s=5)  
axs[0, 0].set_title(f'Eigenvalue {eigenvalue[0]}')
axs[0, 0].set_xlabel('x')
axs[0, 0].set_ylabel('u1')
axs[0, 0].grid(True)

axs[0, 1].scatter(x, u2,s=5)  
axs[0, 1].set_title(f'Eigenvalue {eigenvalue[1]}')
axs[0, 1].set_xlabel('x')
axs[0, 1].set_ylabel('u2')
axs[0, 1].grid(True)

axs[0, 2].scatter(x, u3,s=5)  
axs[0, 2].set_title(f'Eigenvalue {eigenvalue[2]}')
axs[0, 2].set_xlabel('x')
axs[0, 2].set_ylabel('u3')
axs[0, 2].grid(True)

axs[0, 3].scatter(x, u4,s=5)  
axs[0, 3].set_title(f'Eigenvalue {eigenvalue[3]}')
axs[0, 3].set_xlabel('x')
axs[0, 3].set_ylabel('u4')
axs[0, 3].grid(True)

axs[1, 0].scatter(x, u5,s=5)  
axs[1, 0].set_title(f'Eigenvalue {eigenvalue[4]}')
axs[1, 0].set_xlabel('x')
axs[1, 0].set_ylabel('u5')
axs[1, 0].grid(True)

axs[1, 1].scatter(x, u6,s=5)  
axs[1, 1].set_title(f'Eigenvalue {eigenvalue[5]}')
axs[1, 1].set_xlabel('x')
axs[1, 1].set_ylabel('u6')
axs[1, 1].grid(True)

axs[1, 2].scatter(x, u7,s=5)  
axs[1, 2].set_title(f'Eigenvalue {eigenvalue[6]}')
axs[1, 2].set_xlabel('x')
axs[1, 2].set_ylabel('u7')
axs[1, 2].grid(True)

axs[1, 3].scatter(x, u8,s=5)  
axs[1, 3].set_title(f'Eigenvalue {eigenvalue[7]}')
axs[1, 3].set_xlabel('x')
axs[1, 3].set_ylabel('u8')
axs[1, 3].grid(True)


# 调整子图之间的间距
plt.tight_layout()
fig.suptitle(f'Question1,SMR,rho1=0.3sin(PI*x),N=299')
# 显示图形
plt.show()
