# !/usr/bin/env python
# !-*- coding:utf-8 -*-
# !@Time   : 2020/11/5 19:13
# !@Author : DongHan Yang
# !@File   : Circuit_BFS3.py

import cirq
from cirq import LineQubit
import random
import numpy as np
import sys
import queue
from sympy import *

sys.setrecursionlimit(10000)

# result_dict结果真值表 key:真值表;value:[代价,node]
result_dict = {}
reserve_dict = {}
# sympy 的符号运算
u11 = Symbol('u11')
u12 = Symbol('u12')
u21 = Symbol('u21')
u22 = Symbol('u22')
U = Matrix([[u11, u12], [u21, u22]])
u0 = np.mat([1, 0])
u1 = np.mat([0, 1])
U00 = np.kron(u0, u0.T)
U11 = np.kron(u1, u1.T)
U01 = np.kron(u0, u1.T)
U10 = np.kron(u1, u0.T)
# U00 = np.mat([[1, 0], [0, 0]])  # 二值逻辑|0><0|
# U11 = np.mat([[0, 0], [0, 1]])  # 二值逻辑|1><1|
# 使用到的门的矩阵
NOT = np.mat([[0, 1], [1, 0]])
Z = np.mat([[1, 0], [0, -1]])
H = (1 / sqrt(2)) * np.mat([[1, 1], [1, -1]])
S = np.mat([[1, 0], [0, I]])
Sdg = np.mat([[1, 0], [0, -I]])
T = np.mat([[1, 0], [0, (1 + I) / sqrt(2)]])
Tdg = np.mat([[1, 0], [0, (1 - I) / sqrt(2)]])  # 矩阵T的共轭转置（Conjugate Transpose）
V = np.dot(np.mat([
    [1, -I],
    [-I, 1]
]), (1 + I) / 2)
# Swap = np.mat([
#     [1, 0, 0, 0],
#     [0, 0, 1, 0],
#     [0, 1, 0, 0],
#     [0, 0, 0, 1]
# ])
Swap = np.kron(U00, U00) + np.kron(U01, U10) + np.kron(U10, U01) + np.kron(U11, U11)
# 查找电路，Toffi,C_V
Toffi = np.mat([[1, 0, 0, 0, 0, 0, 0, 0]
                   , [0, 1, 0, 0, 0, 0, 0, 0]
                   , [0, 0, 1, 0, 0, 0, 0, 0]
                   , [0, 0, 0, 1, 0, 0, 0, 0]
                   , [0, 0, 0, 0, 1, 0, 0, 0]
                   , [0, 0, 0, 0, 0, 1, 0, 0]
                   , [0, 0, 0, 0, 0, 0, 0, 1]
                   , [0, 0, 0, 0, 0, 0, 1, 0]])
Result_unitary = Swap


# 对角矩阵
def e(n):
    return np.mat(np.diag([1] * (2 ** n)))


# 控制V门
# [[1,0 ,0, 0],
#  [0, 1, 0, 0],
#  [0, 0, 1/2 + I/2, -I*(1/2 + I/2)],
#  [0, 0, -I*(1/2 + I/2), 1/2 + I/2]]
C_V = np.kron(U00, e(1)) + np.kron(U11, V)


# 单量子门算矩阵
def matrix_type(i, n, matrix1):
    for k in range(n):
        if k < i:
            matrix1 = np.kron(e(1), matrix1)
        if k > i:
            matrix1 = np.kron(matrix1, e(1))
    return matrix1


# 符号化简
def simply_unity(cir):
    cir_array = cir.getA()  # 将一个矩阵转化为数组
    # 遍历数组进行化简
    for i in range(len(cir_array)):
        for j in range(len(cir_array[i])):
            cir_array[i][j] = simplify(cir_array[i][j])
    # print(cir_array)
    return cir_array


# ---@----      |0><0|@I + |1><1|@ X
# ---|----
# ---X----
# j受控点,i控制点;计算矩阵
# ---X----      I@|0><0| + X@|1><1|
# ---|----
# ---@----
def matrix_type1(i, j, n, matrix1) -> Matrix:
    cal_num = 0
    for k in range(n):
        k = n - k - 1
        if k != i and k != j and k < j:
            cal_num += 1
            matrix1 = np.kron(e(1), matrix1)
        elif k != i and k != j and k > j:
            cal_num += 1
            matrix1 = np.kron(matrix1, e(1))
        elif k == i and i < j:
            cal_num += 1
            matrix1 = np.kron(U00, e(cal_num)) + np.kron(U11, matrix1)
        elif k == i and i > j:
            cal_num += 1
            matrix1 = np.kron(e(cal_num), U00) + np.kron(matrix1, U11)
    return matrix1


def matrix_type2(i, j, n, matrix1) -> Matrix:
    time = 2
    matrix2 = e(n)
    matrix3 = matrix1
    while time:
        cal_num = 0
        matrix1 = matrix3
        for k in range(n):
            k = n - k - 1
            if k != i and k != j and k < j:
                cal_num += 1
                matrix1 = np.kron(e(1), matrix1)
            elif k != i and k != j and k > j:
                cal_num += 1
                matrix1 = np.kron(matrix1, e(1))
            elif k == i and i < j:
                cal_num += 1
                matrix1 = np.kron(U00, e(cal_num)) + np.kron(U11, matrix1)
            elif k == i and i > j:
                cal_num += 1
                matrix1 = np.kron(e(cal_num), U00) + np.kron(matrix1, U11)
        time -= 1
        i, j = j, i
        matrix2 = np.dot(matrix1, matrix2)
    # print(np.dot(matrix2, matrix1))
    return np.dot(matrix2, matrix1)


# father_node:父节点;children_list子节点;my_type:元电路对象;unitary矩阵,child=0减枝
class EveryNode:
    def __init__(self, father_node, my_type, cost, unitary):
        self.father_node = father_node
        self.children_list = []
        self.my_type = my_type
        self.cost = cost
        self.unitary = unitary
        self.child = 1

    # 添加父亲的子节点,目的为了减枝
    def add_children(self, children):
        self.children_list.append(children)

    # 当发现更加优秀的子孙,老的需要斩掉所有子孙,child置0,列表清空
    def clear_children(self):
        self.child = 0
        for children in self.children_list:
            children.clear_children()
            children.children_list.clear()


# 矩阵-->数字列表(真值表存储)(未使用)
def unitary_list(unitary):
    code_list = [j for unitary_i in unitary for j in range(len(unitary_i)) if unitary_i[j] == 1]
    return code_list


# 元电路生成[电路,代价,矩阵];X,CX,CCX
def type_circuit(n):
    type_x = [[cirq.X(LineQubit(i)), 1, matrix_type(i, n, NOT)] for i in range(n)]
    type_h = [[cirq.H(LineQubit(i)), 1, matrix_type(i, n, H)] for i in range(n)]
    type_z = [[cirq.Z(LineQubit(i)), 1, matrix_type(i, n, Z)] for i in range(n)]
    type_s = [[cirq.S(LineQubit(i)), 1, matrix_type(i, n, S)] for i in range(n)]
    type_sdg = [[(cirq.S ** -1)(LineQubit(i)), 1, matrix_type(i, n, Sdg)] for i in range(n)]
    type_t = [[cirq.T(LineQubit(i)), 1, matrix_type(i, n, T)] for i in range(n)]
    type_tdg = [[(cirq.T ** -1)(LineQubit(i)), 1, matrix_type(i, n, Tdg)] for i in range(1, n)]
    type_cx = [[cirq.CNOT(LineQubit(i), LineQubit(j)), 3, matrix_type1(i, j, n, NOT)] for i in range(n) for j in
               range(n) if i != j]
    type_cz = [[cirq.CZ(LineQubit(i), LineQubit(j)), 3, matrix_type1(i, j, n, Z)] for i in range(n) for j in
               range(n) if i != j]
    type_swap = [[cirq.SWAP(LineQubit(i), LineQubit(j)), 9, matrix_type2(i, j, n, NOT)] for i in range(n) for j in
                 range(i+1, n)]
    type_all = type_x + type_h + type_z + type_s + type_t + type_sdg + type_tdg + type_cx + type_cz + type_swap
    # type_all = type_cx
    print(type_all)
    show_circuit(type_all)
    print('一共{}种情况\n'.format(len(type_all)))
    return type_all


# 显示原电路
def show_circuit(type_all):
    for i, circuit_i in enumerate(type_all):
        circuit = cirq.Circuit(circuit_i[0])
        print("第", i, "情况:\n")
        print(circuit, "\n")


# 根据结果,在结果里面查找,生成电路
def make_circuit(truth_tuple, type_all, result_dict1):
    circuit = []
    min_total_cost = 0
    if truth_tuple in result_dict1:
        node_re = result_dict1[truth_tuple][1]
        min_total_cost = result_dict1[truth_tuple][0]
        while node_re is not None and node_re.my_type is not None:
            type_index = node_re.my_type
            my_type = type_all[type_index][0]
            circuit.append(my_type)
            node_re = node_re.father_node
        # if (node_re.unitary != Swap_unitary).any():
        #     circuit
    return circuit[::-1], min_total_cost


def unitary_to_tuple(circuit_new_unitary):
    unitary_tuple_new = [i for item in circuit_new_unitary for i in item]
    unitary_tuple_new = tuple(unitary_tuple_new)
    return unitary_tuple_new


# 广度优先,循环寻找最优解(递归)
# type_all:基本元情况； node:父节点
def com_circuit(plies_node, plies_node1, type_all):
    while not plies_node.empty() or not plies_node1.empty():
        plies_node, flag = two_queue(plies_node, type_all, result_dict, reserve_dict)
        if flag:
            return
        plies_node1, flag = two_queue(plies_node1, type_all, reserve_dict, result_dict)
        if flag:
            return


def two_queue(plies_node, type_all, result_dict1, reserve_dict1):
    plies_node_list_size = plies_node.qsize()
    # 双向广度优先，每层遍历
    for i in range(plies_node_list_size):
        node = plies_node.get()
        if node.child == 0:  # 节点非最优
            continue
        for circuit_i in range(len(type_all)):
            if circuit_i == node.my_type:  # 相同的会抵消,不考虑
                continue
            else:
                # 优化后 自己算的矩阵
                circuit_new_unitary = simply_unity(np.dot(type_all[circuit_i][2], node.unitary))
                unitary_tuple_new = unitary_to_tuple(circuit_new_unitary)
                node_cost = node.cost + type_all[circuit_i][1]
                node_new = EveryNode(node, circuit_i, node_cost, circuit_new_unitary)
                node.add_children(node_new)
                if unitary_tuple_new in reserve_dict1:
                    plies_node.put(node_new)
                    result_dict1[unitary_tuple_new] = [node_cost, node_new]
                    circuit = cirq.Circuit()
                    # 修改非result的查找
                    circuit_list, min_total_cost = make_circuit(unitary_tuple_new, type_all, result_dict1)
                    circuit_list1, min_total_cost1 = make_circuit(unitary_tuple_new, type_all, reserve_dict1)
                    print(circuit_list)
                    print(circuit_list1)
                    circuit_list += circuit_list1
                    min_total_cost += min_total_cost1
                    circuit.append(circuit_list)
                    check_result_unitary = circuit.unitary()
                    check_result = (check_result_unitary == Result_unitary).all()
                    print(circuit)
                    print('此真正表最低代价为：{}'.format(min_total_cost))
                    print('验证电路矩阵为：\n{}'.format(check_result_unitary))
                    print('初始电路矩阵为：\n{}'.format(Result_unitary))
                    print('验证矩阵相等吗:{}'.format(check_result))
                    return plies_node, 1
                if unitary_tuple_new not in result_dict1:  # 新的解答
                    # print('目前解数目：{}'.format(len(result_dict)))
                    plies_node.put(node_new)
                    result_dict1[unitary_tuple_new] = [node_cost, node_new]
                elif unitary_tuple_new in result_dict1 and result_dict1[unitary_tuple_new][0] > node_cost:  # 更优解
                    plies_node.put(node_new)
                    node_get = result_dict1[unitary_tuple_new][1]  # 找到非优解进行child=0,child_list置空操作
                    node_get.clear_children()
                    result_dict1[unitary_tuple_new] = [node_cost, node_new]  # 更新真值字典
                else:
                    node_new.clear_children()
                    continue
    return plies_node, 0


# 初始化
def init_truth(n):
    circuit_init_unitary = np.identity(2 ** n)
    node_0 = EveryNode(None, None, 0, circuit_init_unitary)
    print(Swap)
    result_tuple = unitary_to_tuple(Result_unitary.getA().tolist())
    print(result_tuple)
    reserve_dict[result_tuple] = [0, None]
    node_1 = EveryNode(None, None, 0, Result_unitary)
    tree_queue = queue.Queue()
    tree_queue.put(node_0)
    tree_queue1 = queue.Queue()
    tree_queue1.put(node_1)
    com_circuit(tree_queue, tree_queue1, type_circuit(n))
    print('(数组列表):[代价，node标号]')
    # print(result_dict)
    print('目前遍历的得到情况数目:{}'.format(len(result_dict)))


def main():
    n = input('请输入电路线路数：\n')
    n = int(n)
    init_truth(n)


if __name__ == '__main__':
    main()
