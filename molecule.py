import math
import random
import time
import openpyxl

# 常量设置 国际单位
sigma = 3.405e-10          # 粒子直径 长度测量尺度
epsilon = 1.654e-21        # 能量尺度
k_B = 1.380649e-23         # 玻尔兹曼常量
mass = 6.69e-26            # 粒子质量
Teq = 400                   # 设定平衡态温度
T0 = Teq/119.8              # 归一化温度
box = [6.0, 3.0, 6.0, 3.0] # 空腔边界 分别为 x 0.5x y 0.5y
N = 20001                  # 执行时间步数，设置为10001是为了后续采样计算平均温度方便
dt = 0.005                  # 时间步长
M_PI = math.pi             # 圆周率
#如果原子距离小于一个sigma，则认为两个原子重叠，设置阈值舍去此次计算
cutoff = 1.1 
cutoffSquare = cutoff * cutoff

# 原子类 代表每个原子的性质
class Atom:
    def __init__(self):
        self.x = [0.0] * N  # 各粒子x坐标
        self.y = [0.0] * N  # 各粒子y坐标
        self.v_x = [0.0] * N  # 各粒子x方向速度
        self.v_y = [0.0] * N  # 各粒子y方向速度
        self.v2 = [0.0] * N   # 各粒子速度平方，方便计算动能
        self.f_x = [0.0] * N  # 各粒子x方向受力
        self.f_y = [0.0] * N  # 各粒子y方向受力
        self.a_x = [0.0] * N  # 各粒子x方向加速度
        self.a_y = [0.0] * N  # 各粒子y方向加速度

# 创建原子数组
atoms = [Atom() for _ in range(16)]
#创建温度数组，用于记录每个时间步长下的温度
T = [0.0] * N
# 设置随机种子
random.seed(time.time())

#最小像约定
def minImageOne(x, box):
    if x > box[0] / 2.0:
        x -= box[0]
    elif x < -box[0] / 2.0:
        x += box[0]
    return x

def minImageTwo(x, y, box):
    x = minImageOne(x, box)
    y = minImageOne(y, box)
    return x, y

#周期性边界处理函数
def pbc(x):
    if x > box[0]:
        x -= box[0]
    elif x < 0:
        x += box[0]
    return x

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
        v0 = 2.9  # 初速度幅值
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

    #计算n=0时每个粒子受力，并更新加速度
    for i in range(16):
        for j in range(i+1, 16):
            xij=atoms[i].x[0]-atoms[j].x[0]
            yij=atoms[i].y[0]-atoms[j].y[0]
            #最小像约定处理
            xij=minImageOne(xij,box)
            yij=minImageOne(yij,box)
            rij=math.sqrt(xij*xij+yij*yij)
            f=24*rij**(-7)-48*rij**(-13)
            atoms[i].f_x[0] += f*xij/rij
            atoms[i].f_y[0] += f*yij/rij
            atoms[j].f_x[0] -= f*xij/rij
            atoms[j].f_y[0] -= f*yij/rij

    #在归一化条件下，粒子加速度和受力相等，所以直接赋值
    for i in range(16):
        atoms[i].a_x[0]=atoms[i].f_x[0]
        atoms[i].a_y[0]=atoms[i].f_y[0]

    #进入n从1到999的循环，每次循环更新速度、位置、加速度
    for n in range(1, N):
        for i in range(16):
            #先更新位置
            atoms[i].x[n] = atoms[i].x[n-1] + atoms[i].v_x[n-1] * dt + 0.5 * atoms[i].a_x[n-1] * dt * dt
            atoms[i].y[n] = atoms[i].y[n-1] + atoms[i].v_y[n-1] * dt + 0.5 * atoms[i].a_y[n-1] * dt * dt

            #再进行周期性边界处理
            atoms[i].x[n] = pbc(atoms[i].x[n])
            atoms[i].y[n] = pbc(atoms[i].y[n])

            #计算受力
            for j in range(i+1, 16):
                xij=atoms[i].x[n]-atoms[j].x[n]
                yij=atoms[i].y[n]-atoms[j].y[n]
                #最小像约定处理
                xij=minImageOne(xij,box)
                yij=minImageOne(yij,box)
                rij=math.sqrt(xij*xij+yij*yij)
                #舍弃粒子重叠的情况
                if rij > cutoff: #如果距离大于cutoff，则计算受力
                    f=24*rij**(-7)-48*rij**(-13)
                    atoms[i].f_x[n] += f*xij/rij
                    atoms[i].f_y[n] += f*yij/rij
                    atoms[j].f_x[n] -= f*xij/rij
                    atoms[j].f_y[n] -= f*yij/rij
                else:#如果距离小于cutoff，则这两个粒子间的作用力为0
                    atoms[i].f_x[n] += 0
                    atoms[i].f_y[n] += 0

            #再更新加速度
            atoms[i].a_x[n] = atoms[i].f_x[n]
            atoms[i].a_y[n] = atoms[i].f_y[n]

            #再更新速度
            atoms[i].v_x[n] = atoms[i].v_x[n-1] + 0.5 * (atoms[i].a_x[n-1] + atoms[i].a_x[n]) * dt
            atoms[i].v_y[n] = atoms[i].v_y[n-1] + 0.5 * (atoms[i].a_y[n-1] + atoms[i].a_y[n]) * dt
            #计算速度平方
            atoms[i].v2[n] = atoms[i].v_x[n]**2 + atoms[i].v_y[n]**2
        
        #计算系统温度，用所有粒子动能的平均值来表示
        T[n] = 0.5*sum(atoms[i].v2[n] for i in range(16))/16

        #判断是否达到平衡，用一段时间内的温度平均值和设定平衡温度的差来判断
        #如果差值比小于0.1，则认为达到平衡，跳出循环
        #一段时间步数为100，如果时间步数是100的倍数，则计算温度平均值
        if n % 500 == 0:
            T_avg = sum(T[n-500:n]) / 500 #计算从0到n的平均值
            if abs(T_avg - T0)/T0 < 0.01:
                #若达到平衡，输出每个粒子的最后一个时间步的速度、温度、对应的步数，写入文件
                #第一列为粒子序号，第二列为步数，第三列为速度，第四列为温度
                workBook = openpyxl.Workbook()
                sheet = workBook.active
                sheet['A1'] = 'atoms'
                sheet['B1'] = 'steps'
                sheet['C1'] = 'v'
                sheet['D1'] = 'T'
                for i in range(16):
                    sheet.append([i, n, math.sqrt(atoms[i].v2[n]), T_avg])
                workBook.save('output.xlsx')
                break
        #如果未达到平衡，则缩放更新速度，继续循环    
            else:
                for i in range(16):
                    atoms[i].v_x[n] = atoms[i].v_x[n] * math.sqrt(T0/T[n])
                    atoms[i].v_y[n] = atoms[i].v_y[n] * math.sqrt(T0/T[n])

def simulateMore():
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
        v0 = 2.9  # 初速度幅值
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

    #计算n=0时每个粒子受力，并更新加速度
    for i in range(16):
        for j in range(i+1, 16):
            xij=atoms[i].x[0]-atoms[j].x[0]
            yij=atoms[i].y[0]-atoms[j].y[0]
            #最小像约定处理
            xij=minImageOne(xij,box)
            yij=minImageOne(yij,box)
            rij=math.sqrt(xij*xij+yij*yij)
            f=24*rij**(-7)-48*rij**(-13)
            atoms[i].f_x[0] += f*xij/rij
            atoms[i].f_y[0] += f*yij/rij
            atoms[j].f_x[0] -= f*xij/rij
            atoms[j].f_y[0] -= f*yij/rij

    #在归一化条件下，粒子加速度和受力相等，所以直接赋值
    for i in range(16):
        atoms[i].a_x[0]=atoms[i].f_x[0]
        atoms[i].a_y[0]=atoms[i].f_y[0]

    #进入n从1到999的循环，每次循环更新速度、位置、加速度
    for n in range(1, N):
        for i in range(16):
            #先更新位置
            atoms[i].x[n] = atoms[i].x[n-1] + atoms[i].v_x[n-1] * dt + 0.5 * atoms[i].a_x[n-1] * dt * dt
            atoms[i].y[n] = atoms[i].y[n-1] + atoms[i].v_y[n-1] * dt + 0.5 * atoms[i].a_y[n-1] * dt * dt

            #再进行周期性边界处理
            atoms[i].x[n] = pbc(atoms[i].x[n])
            atoms[i].y[n] = pbc(atoms[i].y[n])

            #计算受力
            for j in range(i+1, 16):
                xij=atoms[i].x[n]-atoms[j].x[n]
                yij=atoms[i].y[n]-atoms[j].y[n]
                #最小像约定处理
                xij=minImageOne(xij,box)
                yij=minImageOne(yij,box)
                rij=math.sqrt(xij*xij+yij*yij)
                #舍弃粒子重叠的情况
                if rij > cutoff: #如果距离大于cutoff，则计算受力
                    f=24*rij**(-7)-48*rij**(-13)
                    atoms[i].f_x[n] += f*xij/rij
                    atoms[i].f_y[n] += f*yij/rij
                    atoms[j].f_x[n] -= f*xij/rij
                    atoms[j].f_y[n] -= f*yij/rij
                else:#如果距离小于cutoff，则这两个粒子间的作用力为0
                    atoms[i].f_x[n] += 0
                    atoms[i].f_y[n] += 0

            #再更新加速度
            atoms[i].a_x[n] = atoms[i].f_x[n]
            atoms[i].a_y[n] = atoms[i].f_y[n]

            #再更新速度
            atoms[i].v_x[n] = atoms[i].v_x[n-1] + 0.5 * (atoms[i].a_x[n-1] + atoms[i].a_x[n]) * dt
            atoms[i].v_y[n] = atoms[i].v_y[n-1] + 0.5 * (atoms[i].a_y[n-1] + atoms[i].a_y[n]) * dt
            #计算速度平方
            atoms[i].v2[n] = atoms[i].v_x[n]**2 + atoms[i].v_y[n]**2
        
        #计算系统温度，用所有粒子动能的平均值来表示
        T[n] = 0.5*sum(atoms[i].v2[n] for i in range(16))/16

        #判断是否达到平衡，用一段时间内的温度平均值和设定平衡温度的差来判断
        #如果差值比小于0.1，则认为达到平衡，跳出循环
        #一段时间步数为100，如果时间步数是100的倍数，则计算温度平均值
        if n % 500 == 0:
            T_avg = sum(T[n-500:n]) / 500 #计算从0到n的平均值
            if abs(T_avg - T0)/T0 < 0.01:
                #若达到平衡，输出每个粒子的最后一个时间步的速度、温度、对应的步数，写入文件
                #第一列为粒子序号，第二列为步数，第三列为速度，第四列为温度
                workBook = openpyxl.load_workbook('output.xlsx')
                sheet = workBook.active
                for i in range(16):
                    sheet.append([i, n, math.sqrt(atoms[i].v2[n]), T_avg])
                workBook.save('output.xlsx')
                break
        #如果未达到平衡，则缩放更新速度，继续循环    
            else:
                for i in range(16):
                    atoms[i].v_x[n] = atoms[i].v_x[n] * math.sqrt(T0/T[n])
                    atoms[i].v_y[n] = atoms[i].v_y[n] * math.sqrt(T0/T[n])

#进行第一次
simulateOne()
#进行10次
for k in range(150):
    #重新归零初始化
    for i in range(16):
        atoms[i].v_x = [0] * N
        atoms[i].v_y = [0] * N
        atoms[i].x = [0] * N
        atoms[i].y = [0] * N
        atoms[i].a_x = [0] * N
        atoms[i].a_y = [0] * N
        atoms[i].f_x = [0] * N
        atoms[i].f_y = [0] * N
    T = [0.0] * N
    simulateMore()
print("运行完毕，数据写入完成！")