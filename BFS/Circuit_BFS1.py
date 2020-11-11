# !/usr/bin/env python
# !-*- coding:utf-8 -*-
# !@Time   : 2020/11/4 20:43
# !@Author : DongHan Yang
# !@File   : Circuit_BFS1.py

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
def make_circuit(truth_tuple, type_all):
    circuit = []
    min_total_cost = 0
    if truth_tuple in result_dict:
        node_re = result_dict[truth_tuple][1]
        min_total_cost = result_dict[truth_tuple][0]
        while node_re is not None and node_re.my_type is not None:
            type_index = node_re.my_type
            my_type = type_all[type_index][0]
            circuit.append(my_type)
            node_re = node_re.father_node
    return circuit[::-1], min_total_cost


# 广度优先,循环寻找最优解
# type_all:基本元情况； node:父节点
def com_circuit(plies_node, type_all):
    while not plies_node.empty():
        plies_node_list_size = plies_node.qsize()
        # 广度优先，每层遍历
        for i in range(plies_node_list_size):
            node = plies_node.get()
            if node.child == 0:  # 节点非最优
                continue
            for circuit_i in range(len(type_all)):
                if circuit_i == node.my_type:  # 相同的会抵消,不考虑
                    continue
                else:
                    # 优化后 自己算的矩阵
                    circuit_new_unitary = np.dot(type_all[circuit_i][2], node.unitary)
                    unitary_tuple_new = tuple(unitary_list(circuit_new_unitary))
                    node_cost = node.cost + type_all[circuit_i][1]
                    node_new = EveryNode(node, circuit_i, node_cost, circuit_new_unitary)
                    node.add_children(node_new)
                    if unitary_tuple_new not in result_dict:  # 新的解答
                        print('目前解数目：{}'.format(len(result_dict)))
                        plies_node.put(node_new)
                        result_dict[unitary_tuple_new] = [node_cost, node_new]
                    elif unitary_tuple_new in result_dict and result_dict[unitary_tuple_new][0] > node_cost:  # 更优解
                        plies_node.put(node_new)
                        node_get = result_dict[unitary_tuple_new][1]  # 找到非优解进行child=0,child_list置空操作
                        node_get.clear_children()
                        result_dict[unitary_tuple_new] = [node_cost, node_new]  # 更新真值字典
                    else:
                        node_new.clear_children()
                        continue


# 初始化
def init_truth(n, type_all):
    circuit_unitary = np.identity(2 ** n)
    node_0 = EveryNode(None, None, 0, circuit_unitary)
    tree_queue = queue.Queue()
    tree_queue.put(node_0)
    com_circuit(tree_queue, type_all)
    print('(数组列表):[代价，node标号]')
    # print(result_dict)
    print('真值表可表达的情况数目:{}'.format(len(result_dict)))


def main():
    n = input('请输入电路线路数：\n')
    n = int(n)
    type_all = type_circuit(n)
    init_truth(n, type_all)
    circuit = cirq.Circuit()
    # 生成随机测试
    my_truth = [i for i in range(2 ** n)]
    random.shuffle(my_truth)
    # 寻找字典里找最优解
    circuit_list, min_total_cost = make_circuit(tuple(my_truth), type_all)
    while len(circuit_list) == 0:
        print('随机生成的输入：{}在目前去情况下无法生成电路'.format(my_truth))
        random.shuffle(my_truth)
        circuit_list, min_total_cost = make_circuit(tuple(my_truth))
    circuit.append(circuit_list)
    # 与cirq的校验
    circuit_unitary = circuit.unitary()
    check_unitary = unitary_list(circuit_unitary)
    print(circuit)
    print(circuit_unitary)
    print('随机生成的输入为：{}'.format(my_truth))
    print('此真正表最低代价为：{}'.format(min_total_cost))
    print('校验真值表：{}'.format(check_unitary))


if __name__ == '__main__':
    main()
