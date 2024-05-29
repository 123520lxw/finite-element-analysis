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
        self.x = x  # x坐标
        self.y = y  # y坐标
        self.no = no  # 节点编号
        self.fixed = False  # 是否被固定


class Element:
    def __init__(self, node):  # 三节点三角形单元
        self.area = 0.0  # 面积
        self.node = node  # 节点集合(包含三个节点)

    def b_array(self):  # 应变矩阵
        nodes = np.array([[node.x, node.y] for node in self.node])  # 将每个节点的x,y坐标提取成一个列表
        x1, y1 = nodes[0]  # 第一个节点x,y坐标
        x2, y2 = nodes[1]  # 第二个节点x,y坐标
        x3, y3 = nodes[2]  # 第三个节点x,y坐标
        self.area = 0.5 * abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))  # 根据三个点坐标求面积
        b = np.array([[y2 - y3, 0, y3 - y1, 0, y1 - y2, 0],
                      [0, x3 - x2, 0, x1 - x3, 0, x2 - x1],
                      [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2]]) / (2 * self.area)
        return b
        # 根据应变矩阵的表达式求出应变矩阵


class Plane:
    def __init__(self, e, nu, l, w, t, x_num, y_num):
        self.displacement = None  # 位移矩阵
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
        self.generate_mesh()  # 划分网格
        self.k_global()  # 求全局刚度矩阵
        self.apply_boundary()  # 施加荷载

    def d_array(self):  # 根据弹性矩阵的表达式求出弹性矩阵
        d = self.e / (1 - self.nu ** 2) * np.array([[1, self.nu, 0],
                                                    [self.nu, 1, 0],
                                                    [0, 0, (1 - self.nu) / 2]])
        return d

    def k_global(self):  # 全局刚度矩阵
        sum_nodes = self.y_num * self.x_num  # 节点数量
        self.k = np.zeros((2*sum_nodes, 2*sum_nodes))  # 先将全局刚度矩阵设为全零矩阵
        for element in self.elements:  # 先求出每个单元的刚度矩阵
            b = element.b_array()  # 应变矩阵
            a = element.area  # 单元面积
            d = self.d_array()  # 弹性矩阵
            k_local = a*np.dot(b.T, np.dot(d, b))*self.t  # 单元刚度矩阵
            for i in range(3):  # 表示单元刚度矩阵的行号
                for j in range(3):  # 表示单元刚度矩阵的列号
                    for ii in range(2):  # 表示单元刚度矩阵的x方向的刚度
                        for jj in range(2):  # 表示单元刚度矩阵的y方向的刚度
                            row = 2 * element.node[i].no + ii  # 全局刚度矩阵上对应的行位置
                            col = 2 * element.node[j].no + jj  # 全局刚度矩阵上对应的列位置
                            self.k[row, col] += k_local[2 * i + ii, 2 * j + jj]  # 将单元刚度矩阵置入全局刚度矩阵

    def generate_mesh(self):
        dx = self.l / (self.x_num-1)  # 三角形单元的底
        dy = self.w / (self.y_num-1)  # 三角形单元的高
        no = 0  # 节点的编号(从左到右,从上到下)
        for j in range(self.y_num):  # 生成节点
            for i in range(self.x_num):
                x = i * dx  # 节点的x坐标
                y = j * dy  # 节点的y坐标
                self.nodes.append(Node(x, y, no))  # 将节点置入平板的节点内
                no += 1  # 编号自增
        for j in range(self.y_num-1):  # 生成单元,将每四个节点分为两个三角形单元
            for i in range(self.x_num-1):
                n1 = i + j * (self.x_num + 1) - j
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

    def apply_load(self, nodes_force):  # 定义荷载,nodes_force为节点,x方向荷载, y方向荷载
        self.f = np.zeros(2*len(self.nodes))  # 荷载矩阵
        for num in nodes_force:
            node, x, y = [i for i in num]
            self.f[2 * node] = x  # 对应节点上x方向上的荷载
            self.f[2 * node + 1] = y  # 对应节点上y方向上的荷载

    def solve(self):
        # 根据 F = K * X 反推出X,即位移矩阵
        self.displacement = np.linalg.solve(self.k, self.f)

    def plot_deformation(self):  # 将求解出的变形可视化
        plt.figure(figsize=(10, 5))  # 窗口大小
        plt.subplot(1, 2, 1)  # 分为两个视图，左边为未变形形状
        plt.title('未变形形状')
        for element in self.elements:
            raw_data = [[node.x, node.y] for node in element.node]
            plt.gca().add_patch(Polygon(raw_data, closed=True, fill=None, edgecolor='blue'))
            # 将平板以三角形单元的形式输出
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)  # 显示网格

        plt.subplot(1, 2, 2)
        plt.title('变形后形状')
        for element in self.elements:
            deployed_data = [[node.x + self.displacement[2*node.no], node.y + self.displacement[2*node.no + 1]] for node in element.node]
            plt.gca().add_patch(Polygon(deployed_data, closed=True, fill=None, edgecolor='red'))
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.tight_layout()  # 自动调整子图间距
        plt.show()  # 绘制出图像

    def save_result(self, name):  # 保存结果
        x = self.displacement[::2]  # x方向上位移
        y = self.displacement[1::2]  # y方向上位移
        x_label = 'x方向上位移'
        y_label = 'y方向上位移'
        df = pd.DataFrame({
            x_label: x,
            y_label: y
        })
        df.to_excel(f'{name}.xlsx', index=False)


def main():  # 主函数plane_ansys.py
    e = 210e9
    nu = 0.3
    length = 1.0
    weight = 0.1
    thickness = 0.001
    num_x = 11
    num_y = 2
    plate = Plane(e, nu, length, weight, thickness, num_x, num_y)
    nodes_force = [[num_x-1, 0, 10000]]  # 荷载施加在平板右上方节点上
    plate.apply_load(nodes_force)
    plate.solve()
    print(plate.f)
    plate.plot_deformation()
    plate.save_result('result')


if __name__ == '__main__':
    main()
