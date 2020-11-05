# 导入需要的包
import numpy as np

np.set_printoptions(threshold=np.inf)
from sympy import *
import time

# sympy 的符号运算
u11 = Symbol('u11')
u12 = Symbol('u12')
u21 = Symbol('u21')
u22 = Symbol('u22')

# 使用到的门的矩阵
NOT = np.mat([[0, 1], [1, 0]])
Z = np.mat([[1, 0], [0, -1]])
H = (1 / sqrt(2)) * np.mat([[1, 1], [1, -1]])
S = np.mat([[1, 0], [0, I]])
T = np.mat([[1, 0], [0, (1 + I) / sqrt(2)]])
Tct = np.mat([[1, 0], [0, (1 - I) / sqrt(2)]])  # 矩阵T的共轭转置（Conjugate Transpose）
U00 = np.mat([[1, 0], [0, 0]])  # 二值逻辑|0><0|
U11 = np.mat([[0, 0], [0, 1]])  # 二值逻辑|1><1|
U = Matrix([[u11, u12], [u21, u22]])


def circuit_matrix(circuit):
    cir_u = e(len(circuit[0]))
    start_time = time.time()
    for i in circuit:
        cir_u = np.dot(gate_matrix(i), cir_u)
        # print(cir_u)
    end_time = time.time()
    print("计算时间：{}".format(end_time - start_time))
    return cir_u


# 生成2**n阶单位矩阵
def e(n):
    return np.mat(np.diag([1] * (2 ** n)))


# 计算单个量子门的矩阵
def gate_matrix(v):
    global matrix1
    if v[0] in ['t1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8']:
        matrix1 = NOT
    elif v[0] == 'H':
        matrix1 = H
    elif v[0] == 'S':
        matrix1 = S
    elif v[0] == 'T':
        matrix1 = T
    elif v[0] == 'Tct':
        matrix1 = Tct
    elif v[0] == 'U':
        matrix1 = U
    elif v[0] == 'Z':
        matrix1 = Z
    n = 0
    for i in range(1, len(v)):
        n = n + 1
        if v[i] == 1:
            # n = n + 1
            matrix1 = np.kron(e(1), matrix1)
        elif v[i] == 2:
            # n = n + 1
            matrix1 = np.kron(U00, e(n)) + np.kron(U11, matrix1)
        elif v[i] == 3:
            # n = n + 1
            matrix1 = np.kron(matrix1, e(1))
        elif v[i] == 4:
            # n = n + 1
            matrix1 = np.kron(e(n), U00) + np.kron(matrix1, U11)
    return matrix1


# 遍历计算的步骤
def traverse(gates):
    steps = []
    for i in range(1, len(gates)):
        steps.append(traverse_gate(gates[0], gates[i]))
    return steps


def traverse_gate(gate_line, gate):
    step = []
    for i in gate_line:
        if i in gate and i == gate[-1]:
            continue
        elif i in gate and i < gate[-1]:
            step.append(2)
        elif i in gate and i > gate[-1]:
            step.append(4)
        elif (i not in gate) and (i < gate[-1]):
            step.append(1)
        elif (i not in gate) and (i > gate[-1]):
            step.append(3)
    step.append(gate[0])
    step.reverse()
    return step


# 从电路文件中获取有用信息
def get_gates(source_file_path):
    # 读取源文件数据
    file = open(source_file_path, "r+")
    lines = []
    gates = []

    # 将电路源文件中每一行内容分别存放到列表lines中
    for item in file:
        lines.append(item.strip("\n"))  # 去掉行末换行符,将每一行添加到列表中
    # print(lines)

    # 截取第一行信息
    gates.append(lines[0][3:].split(","))

    # 截取 BEGIN 和 END 之间的电路门信息
    for i in range(lines.index("BEGIN") + 1, lines.index("END")):
        gates.append(lines[i].replace(",", " ").split(" "))
        # print(gates)
    return traverse(gates)


# 主函数
def cal_circuit_matrix(source_file_path):
    circuit = get_gates(source_file_path)
    cir = circuit_matrix(circuit)
    cir_array = cir.getA()  # 将一个矩阵转化为数组
    # 遍历数组进行化简
    for i in range(len(cir_array)):
        for j in range(len(cir_array[i])):
            cir_array[i][j] = simplify(cir_array[i][j])
    print(cir_array)


# 程序入口
if __name__ == '__main__':
    # 电路源文件的地址
    circuit_source_file_path = ".\Toffoli.tfc"
    cal_circuit_matrix(circuit_source_file_path)
