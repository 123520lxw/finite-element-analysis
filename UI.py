# -*- coding:utf-8 -*-#
# @Time :2024/5/19 10:49
# @AUTHOR :李雪巍
# @FILE :test.py
# @Software :PyCharm
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import *
from plane_ansys import Plane

class UI:
    def __init__(self, root):
        self.root = root
        self.root.title("平面应力问题")

        ttk.Label(root, text="长度 (m):").grid(row=1, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="宽度 (m):").grid(row=3, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="厚度 (m):").grid(row=5, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="弹性模量 (Pa):").grid(row=7, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="泊松比:").grid(row=9, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="x方向上节点数:").grid(row=11, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="y方向上节点数:").grid(row=13, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="受力节点(以,分隔):").grid(row=15, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="节点对应的力(以;分隔):").grid(row=17, column=0, padx=5, pady=5, sticky="e")
        ttk.Label(root, text="保存文件名:").grid(row=19, column=0, padx=5, pady=5, sticky="e")
        # 输入框
        self.length_entry = ttk.Entry(root)
        self.length_entry.grid(row=1, column=1, padx=5, pady=5, sticky="w")

        self.width_entry = ttk.Entry(root)
        self.width_entry.grid(row=3, column=1, padx=5, pady=5, sticky="w")

        self.thickness_entry = ttk.Entry(root)
        self.thickness_entry.grid(row=5, column=1, padx=5, pady=5, sticky="w")

        self.elastic_modulus_entry = ttk.Entry(root)
        self.elastic_modulus_entry.grid(row=7, column=1, padx=5, pady=5, sticky="w")

        self.poissons_ratio_entry = ttk.Entry(root)
        self.poissons_ratio_entry.grid(row=9, column=1, padx=5, pady=5, sticky="w")

        self.x_num_entry = ttk.Entry(root)
        self.x_num_entry.grid(row=11, column=1, padx=5, pady=5, sticky="w")

        self.y_num_entry = ttk.Entry(root)
        self.y_num_entry.grid(row=13, column=1, padx=5, pady=5, sticky="w")

        self.load_nodes = ttk.Entry(root)
        self.load_nodes.grid(row=15, column=1, padx=5, pady=5, sticky="w")

        self.loads = ttk.Entry(root)
        self.loads.grid(row=17, column=1, padx=5, pady=5, sticky="w")

        self.name = ttk.Entry(root)
        self.name.grid(row=19, column=1, padx=5, pady=5, sticky="w")

        # 提交按钮
        ttk.Button(root, text="确定", command=self.submit).grid(row=20, columnspan=2, padx=5, pady=10)

    def submit(self):
        # 获取输入值
        length = float(self.length_entry.get())
        width = float(self.width_entry.get())
        thickness = float(self.thickness_entry.get())
        elastic_modulus = float(self.elastic_modulus_entry.get())
        poissons_ratio = float(self.poissons_ratio_entry.get())
        x_num = int(self.x_num_entry.get())
        y_num = int(self.y_num_entry.get())
        lodes_nodes = self.load_nodes.get().split(',')
        nodes = list(map(int, lodes_nodes))
        loads = self.loads.get().split(';')
        nodes_force = []
        for i in range(len(loads)):
            x, y = list(map(int,loads[i].split(',')))
            node = nodes[i]
            nodes_force.append([node, x, y])
        name = self.name.get()
        # 将输入值传递给主函数
        main(length, width, thickness, elastic_modulus, poissons_ratio, x_num, y_num, nodes_force, name)


def main(length, width, thickness, elastic_modulus, poissons_ratio, x_num, y_num, load, name):
    plate = Plane(elastic_modulus, poissons_ratio, length, width, thickness, x_num, y_num)
    plate.apply_load(load)
    plate.solve()
    fileDir = askdirectory()
    filename = fileDir + '/' + name
    print(filename)
    plate.save_result(filename)
    return


if __name__ == "__main__":
    root = tk.Tk()
    app = UI(root)
    root.mainloop()