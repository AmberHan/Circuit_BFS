# !/usr/bin/env python
# !-*- coding:utf-8 -*-
# !@Time   : 2020/11/4 20:43
# !@Author : DongHan Yang
# !@File   : Circuit_BFS1.1.py

import cirq
from cirq import LineQubit
import random
import numpy as np
import sys
import queue

sys.setrecursionlimit(100000)

# result_dict结果真值表 key:真值表;value:[代价,node]
result_dict = {}


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


# 矩阵-->数字列表(真值表存储)
def unitary_list(unitary):
    code_list = [j for unitary_i in unitary for j in range(len(unitary_i)) if unitary_i[j] == 1]
    return code_list


# 元电路生成[电路,代价,矩阵];X,CX,CCX
def type_circuit(n):
    type_x = [[cirq.X(LineQubit(i)), 1] for i in range(n)]
    type_cnot = [[cirq.CNOT(LineQubit(i), LineQubit(j)), 3] for i in range(n) for j in range(n) if i != j]
    type_ccnot = [[cirq.CCNOT(LineQubit(i), LineQubit(j), LineQubit(k)), 7] for i in range(n)
                  for j in range(n) for k in range(j, n) if i != j and i != k and j != k]
    type_all = type_x + type_cnot + type_ccnot
    # print('{}\n{}种情况'.format(type_all, len(type_all)))
    print('{}种情况'.format(len(type_all)))
    q = [LineQubit(i) for i in range(n)]
    for i in range(len(type_all)):
        circuit = cirq.Circuit(cirq.I.on_each(*q))
        circuit.append(type_all[i][0])
        type_all[i].append(circuit.unitary())
    return type_all


# 根据结果,在结果里面查找,生成电路
def make_circuit(truth_tuple):
    circuit = []
    min_total_cost = 0
    if truth_tuple in result_dict:
        node_re = result_dict[truth_tuple][1]
        min_total_cost = result_dict[truth_tuple][0]
        while node_re is not None:
            my_type = node_re.my_type
            circuit.append(my_type)
            node_re = node_re.father_node
    return circuit[::-1], min_total_cost


# 广度优先,循环寻找最优解
# type_all:基本元情况； node:父节点
def com_circuit(plies_node, type_all):
    while plies_node:
        # 广度优先，每层遍历
        next_plies_node = []
        for node in plies_node:
            if node.child == 0:  # 节点非最优
                continue
            for circuit_i in type_all:
                if circuit_i[0] == node.my_type:  # 相同的会抵消,不考虑
                    continue
                else:
                    # 优化后 自己算的矩阵
                    circuit_new_unitary = np.dot(circuit_i[2], node.unitary)
                    unitary_tuple_new = tuple(unitary_list(circuit_new_unitary))
                    node_cost = node.cost + circuit_i[1]
                    node_new = EveryNode(node, circuit_i[0], node_cost, circuit_new_unitary)
                    node.add_children(node_new)
                    if unitary_tuple_new not in result_dict:  # 新的解答
                        print('目前解数目：{}'.format(len(result_dict)))
                        next_plies_node.append(node_new)
                        result_dict[unitary_tuple_new] = [node_cost, node_new]
                    elif unitary_tuple_new in result_dict and result_dict[unitary_tuple_new][0] > node_cost:  # 更优解
                        next_plies_node.append(node_new)
                        node_get = result_dict[unitary_tuple_new][1]  # 找到非优解进行child=0,child_list置空操作
                        node_get.clear_children()
                        result_dict[unitary_tuple_new] = [node_cost, node_new]  # 更新真值字典
                    else:
                        node_new.clear_children()
                        node_new.child = 0
        plies_node = next_plies_node.copy()


# 初始化
def init_truth(n):
    q = [LineQubit(i) for i in range(n)]
    init_type = cirq.I.on_each(*q)
    circuit_unitary = np.identity(2 ** n)
    node_0 = EveryNode(None, init_type, 0, circuit_unitary)
    # tree_queue = queue.Queue()
    tree_queue = [node_0]
    com_circuit(tree_queue, type_circuit(n))
    print('(数组列表):[代价，node标号]')
    # print(result_dict)
    print('真值表可表达的情况数目:{}'.format(len(result_dict)))


final_min_cost = 30
final_list = []
final_circuit_list = []


def print_circuit(my_truth):
    circuit = cirq.Circuit()
    circuit_list, min_total_cost = make_circuit(tuple(my_truth))
    circuit.append(circuit_list)
    # 与cirq的校验
    circuit_unitary = circuit.unitary()
    check_unitary = unitary_list(circuit_unitary)
    print(circuit)
    # print(circuit_unitary)
    print('随机生成的输入为：{}'.format(my_truth))
    print('此真正表最低代价为：{}'.format(min_total_cost))
    print('校验真值表：{}'.format(check_unitary))


def print_circuits(circuits_list):
    circuit = cirq.Circuit()
    circuit.append(circuits_list)
    # 与cirq的校验
    circuit_unitary = circuit.unitary()
    check_unitary = unitary_list(circuit_unitary)
    print(circuit)
    # print(circuit_unitary)
    print('随机生成的输入为：{}'.format(final_list))
    print('此真正表最低代价为：{}'.format(final_min_cost))
    print('校验真值表：{}'.format(check_unitary))


def find_circuit(my_truth):
    # 寻找字典里找最优解
    global final_min_cost
    global final_list
    global final_circuit_list
    circuit_list, min_total_cost = make_circuit(tuple(my_truth))
    if min_total_cost < final_min_cost:
        final_min_cost = min_total_cost
        final_circuit_list = circuit_list
        final_list = my_truth
    # while len(circuit_list) == 0:
    #     print('随机生成的输入：{}在目前去情况下无法生成电路'.format(my_truth))
    #     random.shuffle(my_truth)
    #     circuit_list, min_total_cost = make_circuit(tuple(my_truth))
    # circuit.append(circuit_list)
    # 与cirq的校验
    # circuit_unitary = circuit.unitary()
    # check_unitary = unitary_list(circuit_unitary)
    # print(circuit)
    # # print(circuit_unitary)
    # print('随机生成的输入为：{}'.format(my_truth))
    # print('此真正表最低代价为：{}'.format(min_total_cost))
    # print('校验真值表：{}'.format(check_unitary))


total_truths = []
even_index = [0, 1, 2, 4]


def solve(rows, is_visit, one_truth):
    if rows == 8:
        total_truths.append(one_truth)
    else:
        # 0 可放置， 1不能放；获取为0的位
        available_positions = 0x0f & (~is_visit)  # 获取可访问的位数 偶数低四位 0000 1111
        leave_positions = 0xf0 & (~is_visit)  # 获取可访问的位数 奇数高四位  1111 0000
        flag, diff = 0, 0
        # leave_positions = leave_positions1
        if rows not in even_index:
            available_positions, leave_positions = leave_positions, available_positions
            flag, diff = 1, 1
        available_positions1 = available_positions
        while available_positions:
            position = available_positions & (-available_positions)  # 取最后一位1
            next_position = (~position) & available_positions1
            this_index = bin(position - 1).count("1")  # 0001 -> 编号0
            if flag:
                this_index -= 4
            one_truth1 = one_truth[:]
            one_truth1.append(2 * this_index + diff)
            available_positions = available_positions & (available_positions - 1)  # 清空最后一位1
            solve(rows + 1, ~(next_position | leave_positions), one_truth1)


def generate(row, one_truth, even, odd):
    if row == 8:
        total_truths.append(one_truth)
    else:
        length = len(odd)
        if row in even_index:
            length = len(even)
        for index in range(length):
            need = odd[:]
            flag = 0
            if row in even_index:
                flag = 1
                need = even[:]
            one_truth1 = one_truth[:]
            one_truth1.append(need.pop(index))
            if flag:
                generate(row + 1, one_truth1, need, odd[:])
            else:
                generate(row + 1, one_truth1, even, need[:])


def main():
    n = input('请输入电路线路数：\n')
    n = int(n)
    init_truth(n)
    # 生成随机测试
    # my_truth = [i for i in range(2 ** n)]
    # random.shuffle(my_truth)
    # [0,1,2,4]->[0,2,4,6]可为偶数
    # solve(0, 0, [])
    generate(0, [], [0, 2, 4, 6], [1, 3, 5, 7])
    print(total_truths)
    print(len(total_truths))
    for my_truth in total_truths:
        find_circuit(my_truth)
    print_circuits(final_circuit_list)
    print_circuit(final_list)


if __name__ == '__main__':
    main()
