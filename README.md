# 简介🌱：
采用广度优先算法,探索出电路构造矩阵，达到构造矩阵作用
1. BFS1和BFS1。1，是使用numpy计算的矩阵，主要针对传统电路(X,CX,CCX门),采用cirq.unity()验证
2. BFS2，元电路更加复杂(但无CV,CZ,Toffi)，是使用符号矩阵自己算的，针对更加复杂的电路(单项链表)
3、BFS3，和BFS2一样，不过使用的双向广度优先；（新增SWAP,Fredkin门）

# 区别💬
### BFS1：
- 主要解决经典电路遍历所有情况，找到所有真值表，元电路简单，都是0,1真值

### BFS1.1:
- 和1一样，但广度优先搜索，未使用队列，而是使用列表，然后整体清空（个人感觉比队列每次移动有所与优化，但是耗时和1相同）
### BFS2：
- 主要是找电路，寻找某个矩阵的电路构造，找到就停止
- 门构造复杂：X,CX,Z,CZ,H,S,Sdg,T,Tdg,Swap,Fredin(控制交换门);（也可构造V,CV,CCNOT矩阵作为元门，暂未使用，因为cirq没有找到cirq.CV等）
- 变化函数：type_circuit(生成门函数),com_circuit(递归函数)
- 代码目前找到SWap门，请输入：2

# cal_matrix文件：
- 根据Toffoli.tfc计算电路矩阵

# 注意🎉
BFS3：行247以及62，查找和校验内容时候，需要修改Swap；因为双向广度优先，结果矩阵往前找，找对应内容
Result_unitary = Swap；check_result = (Swap == Result_unitary).all()

# 还未修改⚡
- SWAP矩阵应该构造的是错的，还是老问题
- BFS3:结果的node附加时候，需要取对应的逆矩阵，即选择电路时候，S->sdg
- 后期需新增CH,CY等更多的门,[cirq门链接](https://cirq.readthedocs.io/en/stable/docs/interop.html?highlight=CY#Gate-conversion-rules)
- cirq.ControlledGate 不会使用，导致构造不了更多的门