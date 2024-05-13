# -*- coding:utf-8 -*-#
# @Time :2024/5/6 9:53
# @AUTHOR :李雪巍
# @FILE :plane_ansys.py
# @Software :PyCharm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pandas as pd
import matplotlib as mpl
mpl.rcParams["font.family"] = "KaiTi"
mpl.rcParams["axes.unicode_minus"]=False # 正常显示负号


class Node:   # 定义节点，由x,y两个坐标表示
    def __init__(self, x, y, no=-1):
        self.x = x
        self.y = y
        self.no = no  # 节点编号
        self.fixed = False  # 是否被固定


class Element:
    def __init__(self, node):  # 三节点三角形单元
        self.area = 0.0
        self.node = node

    def b_array(self):  # 应变矩阵
        nodes = np.array([[node.x, node.y] for node in self.node])
        x1, y1 = nodes[0]
        x2, y2 = nodes[1]
        x3, y3 = nodes[2]
        self.area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        b = np.array([[y2 - y3, 0, y3 - y1, 0, y1 - y2, 0],
                      [0, x3 - x2, 0, x1 - x3, 0, x2 - x1],
                      [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2]]) / (2 * self.area)
        return b


class Plane:
    def __init__(self, e, nu, l, w, t, x_num, y_num):
        self.displacement = None
        self.elements = []  # 单元集合
        self.nodes = []  # 节点集合
        self.e = e  # 弹性模量
        self.nu = nu  # 泊松比
        self.l = l  # 平板长度
        self.w = w  # 平板宽度
        self.t = t  # 平板厚度
        self.x_num = x_num  # x方向节点数量
        self.y_num = y_num  # y方向节点数量
        self.k = None  # 全局刚度矩阵
        self.f = None  # 荷载矩阵
        self.generate_mesh()
        self.k_global()
        self.apply_boundary()

    def d_array(self):  # 弹性矩阵
        d = self.e / (1 - self.nu ** 2) * np.array([[1, self.nu, 0],
                                                    [self.nu, 1, 0],
                                                    [0, 0, (1 - self.nu) / 2]])
        return d

    def k_global(self):  # 全局刚度矩阵
        sum_nodes = self.y_num * self.x_num
        self.k = np.zeros((2*sum_nodes, 2*sum_nodes))
        for element in self.elements:
            b = element.b_array()
            a = element.area
            d = self.d_array()
            k_local = a*np.dot(b.T, np.dot(d, b))*self.t  # 单元刚度矩阵
            for i in range(3):  # 每个单元有三个节点
                for j in range(3):
                    for ii in range(2):  # 每个节点两个自由度
                        for jj in range(2):
                            row = 2 * element.node[i].no + ii
                            col = 2 * element.node[j].no + jj
                            self.k[row, col] += k_local[2 * i + ii, 2 * j + jj]

    def generate_mesh(self):
        dx = self.l / (self.x_num-1)
        dy = self.w / (self.y_num-1)
        no = 0
        for j in range(self.y_num):  # 生成节点
            for i in range(self.x_num):
                x = i * dx
                y = j * dy
                self.nodes.append(Node(x, y, no))
                no += 1
        for j in range(self.y_num-1):  # 生成单元
            for i in range(self.x_num-1):
                n1 = i + j * (self.x_num + 1)
                n2 = n1 + 1
                n3 = n1 + self.x_num
                n4 = n3 + 1
                nodes1 = [self.nodes[n1], self.nodes[n2], self.nodes[n3]]
                self.elements.append(Element(nodes1))
                nodes2 = [self.nodes[n2], self.nodes[n4], self.nodes[n3]]
                self.elements.append(Element(nodes2))

    def apply_boundary(self):
        for i in range(self.y_num):
            self.k[2 * i * self.x_num, :] = 0  # 节点的x所在的行和列
            self.k[:, 2 * i * self.x_num] = 0
            self.k[2 * i * self.x_num + 1, :] = 0  # 节点的y所在的行和列
            self.k[:, 2 * i * self.x_num + 1] = 0
            self.k[2 * i * self.x_num, 2 * i * self.x_num] = 1
            self.k[2 * i * self.x_num + 1, 2 * i * self.x_num + 1] = 1
            # 将左侧固定的节点的刚度矩阵行和列置零，并设置对角线上的值为1，表示这些自由度被完全约束

    def apply_load(self, nodes, force, direction=1):  # 定义荷载，1表示y方向,0表示x方向
        self.f = np.zeros(2*len(self.nodes))
        for node in nodes:
            self.f[2*node + direction] = force/len(nodes)

    def solve(self):
        self.displacement = np.linalg.solve(self.k, self.f)

    def plot_deformation(self):
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.title('未变形形状')
        for element in self.elements:
            raw_data = [[node.x, node.y] for node in element.node]
            plt.gca().add_patch(Polygon(raw_data, closed=True, fill=None, edgecolor='blue'))
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)

        plt.subplot(1, 2, 2)
        plt.title('变形后形状')
        for element in self.elements:
            deployed_data = [[node.x + self.displacement[2*node.no], node.y + self.displacement[2*node.no + 1]] for node in element.node]
            plt.gca().add_patch(Polygon(deployed_data, closed=True, fill=None, edgecolor='red'))
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)

        plt.tight_layout()
        plt.savefig('./result.jpg')

def load_ansys_data():
    x_displacement = pd.read_excel('./ansys_result.xlsx', sheet_name='Sheet1', header=None)[1]
    y_displacement = pd.read_excel('./ansys_result.xlsx', sheet_name='Sheet2', header=None)[1]
    return x_displacement, y_displacement


e = 210e9
nu = 0.3
length = 1.0
weight = 0.1
thickness = 0.001
num_x = 11
num_y = 2
plate = Plane(e, nu, length, weight, thickness, num_x, num_y)
plate.apply_boundary()
nodes = [num_x-1]  # 荷载施加在平板右上方节点上
plate.apply_load(nodes, 10000)
plate.solve()
plate.plot_deformation()
x1 = plate.displacement[::2]
y1 = plate.displacement[1::2]
data = {
    'x方向上位移': x1,
    'y方向上唯一': y1
}
pd1 = pd.DataFrame(data)
pd1.to_excel('./result.xlsx', index=False)