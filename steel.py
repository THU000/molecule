#本程序为模拟16个原子在边长为6方形空腔中的钢分子，并计算达到平衡态时的平均分子自由程
import math
import random
import time
import numpy as np

# 常量设置 国际单位
sigma = 3.405e-10          # 粒子直径 长度测量尺度
epsilon = 1.654e-21        # 能量尺度
k_B = 1.380649e-23         # 玻尔兹曼常量
mass = 6.69e-26            # 粒子质量
Teq = 300                   # 设定平衡态温度
T0 = Teq/119.8              # 归一化温度
box = [6.0, 3.0, 6.0, 3.0] # 空腔边界 分别为 x 0.5x y 0.5y
N = 40001                  # 执行时间步数，设置为10001是为了后续采样计算平均温度方便
dt = 0.005                  # 时间步长
M_PI = math.pi             # 圆周率

# 原子类 代表每个原子的性质
class Atom:
    def __init__(self):
        self.x = [0.0] * N  # 各粒子x坐标
        self.y = [0.0] * N  # 各粒子y坐标
        self.v_x = [0.0] * N  # 各粒子x方向速度
        self.v_y = [0.0] * N  # 各粒子y方向速度
        self.v2 = [0.0] * N   # 各粒子速度平方，方便计算动能
        self.s = [0.0] * N  # 各粒子在每一个时间步的位移
        self.collisions = [0.0]*N          #各粒子经历的碰撞次数

# 创建原子数组
atoms = [Atom() for _ in range(16)]

#创建温度数组，用于记录每个时间步长下的温度
T = [0.0] * N
# 设置随机种子
random.seed(time.time())

#周期性边界处理函数
def pbc(x):
    if x > box[0]:
        x -= box[0]
    elif x < 0:
        x += box[0]
    return x

#交换碰撞方向速度，垂直碰撞方向的速度不变
def calculate_new_velocity(v1, v2,  center1, center2):
    """
    计算碰撞后的新速度。
    
    :param v1: 第一个球的速度，格式为(vx1, vy1)
    :param v2: 第二个球的速度，格式为(vx2, vy2)
    :param collision_point: 碰撞点的坐标
    :param center1: 第一个球的球心坐标
    :param center2: 第二个球的球心坐标
    :return: 碰撞后的新速度，格式为((new_vx1, new_vy1), (new_vx2, new_vy2))
    """
    # 计算球心连线方向的单位向量
    direction_vector = (center2[0] - center1[0], center2[1] - center1[1])
    distance = math.sqrt(direction_vector[0]**2 + direction_vector[1]**2)
    unit_vector = (direction_vector[0] / distance, direction_vector[1] / distance)
    
    # 计算两个球在球心连线方向上的速度分量
    v1_parallel = dot_product(v1, unit_vector)
    v2_parallel = dot_product(v2, unit_vector)
    
    # 计算垂直于球心连线方向的速度分量
    v1_perpendicular = (v1[0] - v1_parallel * unit_vector[0], v1[1] - v1_parallel * unit_vector[1])
    v2_perpendicular = (v2[0] - v2_parallel * unit_vector[0], v2[1] - v2_parallel * unit_vector[1])
    
    # 碰撞后交换球心连线方向的速度分量
    new_v1_parallel = v2_parallel
    new_v2_parallel = v1_parallel
    
    # 计算新的速度
    new_v1 = (new_v1_parallel * unit_vector[0] + v1_perpendicular[0], new_v1_parallel * unit_vector[1] + v1_perpendicular[1])
    new_v2 = (new_v2_parallel * unit_vector[0] + v2_perpendicular[0], new_v2_parallel * unit_vector[1] + v2_perpendicular[1])
    
    return new_v1, new_v2

def dot_product(v1, v2):
    """
    计算两个向量的点积。
    
    :param v1: 向量1
    :param v2: 向量2
    :return: 点积
    """
    return v1[0] * v2[0] + v1[1] * v2[1]

#单次模拟
def simulateOne():
    # 初始化每个原子的位置
    dx = box[0] / 4.0  # 每个小格子的x长度
    dy = box[2] / 4.0  # 每个小格子的y长度

    for i in range(16):
        row = i // 4  # 计算行号
        col = i % 4   # 计算列号
        atoms[i].x[0] = dx * (col + 0.5)  # 设置x坐标为中心
        atoms[i].y[0] = dy * (row + 0.5)  # 设置y坐标为中心

    # 初始化每个原子的速度，初速度幅值确定为2.9，方向随机
    for i in range(16):
        theta = random.uniform(0, 2 * M_PI)  # 随机角度
        v0 = 2.2  # 初速度幅值
        atoms[i].v_x[0] = v0 * math.cos(theta)
        atoms[i].v_y[0] = v0 * math.sin(theta)
    
    #使质心速度为零
    totalMass = 0.0
    centerOfMassVelocity = [0.0, 0.0]
    for i in range(16):
        totalMass += 1.0
        centerOfMassVelocity[0] += atoms[i].v_x[0]*1.0
        centerOfMassVelocity[1] += atoms[i].v_y[0]*1.0
    centerOfMassVelocity[0] /= totalMass
    centerOfMassVelocity[1] /= totalMass
    for i in range(16):
        atoms[i].v_x[0] -= centerOfMassVelocity[0]
        atoms[i].v_y[0] -= centerOfMassVelocity[1]
        #计算速度平方
        atoms[i].v2[0] = atoms[i].v_x[0]**2 + atoms[i].v_y[0]**2

    #进入n从1到999的循环，每次循环判断是否发生碰撞，如果发生碰撞，则交换碰撞方向的速度，垂直于碰撞方向的速度不变
    #判断粒子是否发生碰撞的条件为粒子距离小于直径
    #碰撞改变速度后，重置
    #对每一个时间步下的每一个粒子都和其他粒子进行碰撞判断，已经判断过的两个粒子不再进行判断
    for n in range(1, N):
        #先进行位置更新
        for i in range(16):
            atoms[i].x[n] = atoms[i].x[n-1] + atoms[i].v_x[n-1] * dt
            atoms[i].y[n] = atoms[i].y[n-1] + atoms[i].v_y[n-1] * dt

            #判断是否越界，越界则周期性边界条件
            atoms[i].x[n] = pbc(atoms[i].x[n])
            atoms[i].y[n] = pbc(atoms[i].y[n])
            #粒子运动距离，速度乘以时间
            atoms[i].s[n] = math.sqrt((atoms[i].v_x[n-1]*dt)**2 + (atoms[i].v_y[n-1]*dt)**2)

        #判断是否发生碰撞
        for i in range(16):
            for j in range(i+1, 16):
                #计算两个粒子之间的距离
                dx = atoms[i].x[n] - atoms[j].x[n]
                dy = atoms[i].y[n] - atoms[j].y[n]

                #判断是否发生碰撞
                if (dx**2 + dy**2) <= 100:
                    v1 = (atoms[i].v_x[n-1], atoms[i].v_y[n-1])
                    v2 = (atoms[j].v_x[n-1], atoms[j].v_y[n-1])
                    center1 = (atoms[i].x[n], atoms[i].y[n])
                    center2 = (atoms[j].x[n], atoms[j].y[n])
                    # 计算碰撞后的新速度
                    new_v1, new_v2 = calculate_new_velocity(v1, v2,  center1, center2)
                    atoms[i].v_x[n] = new_v1[0]
                    atoms[i].v_y[n] = new_v1[1]
                    atoms[j].v_x[n] = new_v2[0]
                    atoms[j].v_y[n] = new_v2[1]

                    #碰撞次数加1
                    atoms[i].collisions[n] = 1.0
                    #print("碰一次")
                else:
                    #速度保持不变
                    atoms[i].v_x[n] = atoms[i].v_x[n-1]
                    atoms[i].v_y[n] = atoms[i].v_y[n-1]

        #计算速度平方
        for i in range(16):
            atoms[i].v2[n] = atoms[i].v_x[n]**2 + atoms[i].v_y[n]**2
        #计算系统温度，用所有粒子动能的平均值来表示
        T[n] = 0.5*sum(atoms[i].v2[n] for i in range(16))/16

        #判断是否达到平衡，用一段时间内的温度平均值和设定平衡温度的差来判断
        #如果差值比小于0.1，则认为达到平衡，跳出循环
        #一段时间步数为500，如果时间步数是500的倍数，则计算温度平均值
        if n % 500 == 0:
            T_avg = sum(T[n-500:n]) / 500 #计算从0到n的平均值
            if abs(T_avg - T0)/T0 < 0.01:
                #总碰撞次数
                k = 0.0
                L = 0.0
                for i in range(16):
                    k += sum(atoms[i].collisions[n-500:n])
                    L += sum(atoms[i].s[n-500:n])
                print("达到平衡，温度为：", T_avg,"\n")
                print(k,"\n")
                print(L,"\n")
                if k == 0:
                    print("没发生碰撞")
                    print(n)
                else:
                    print("平均自由程：", L/k,"\n")
                    print(n)
                break
            else:
                for i in range(16):
                    atoms[i].v_x[n] = atoms[i].v_x[n] * math.sqrt(T0/T[n])
                    atoms[i].v_y[n] = atoms[i].v_y[n] * math.sqrt(T0/T[n])

simulateOne()
print("模拟结束")