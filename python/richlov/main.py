import numpy as np
import matplotlib.pyplot as plt

from z1 import u_xt
from z2 import w_xt


if __name__ == "__main__":
    N = 150  # максимальное число шагов по х
    K = 500  # максимальное число шагов по t
    l = 1  # значение х на правой границе
    h = l / N  # шаг сетки по х
    T = 1  # максимальное значение времени t на правой границе
    t = T / K  # шаг сетки по времени

    # сетка
    x_i = np.arange(0, 1, h)  # значения в узлах по х
    t_j = np.arange(0, 1, t)  # значение в узлах по t
    r_j = len(t_j)  # количество узлов по t
    r_i = len(x_i)  # количество узлов по x
    w_h_t = np.zeros([r_j, r_i])  # итоговая сетка размером x_i*t_j
    a = np.array([1, 0])

    omega1, omega2 = -1, 2

    res1 = u_xt(w_h_t.copy(), x_i, t_j, (omega1, omega2))
    res2 = w_xt(w_h_t.copy(), x_i, t_j, (omega1, omega2))
    diff = np.abs(res2) - np.abs(res1)
    print(f"Максимальная разность элементов в массивах = {np.max(np.abs(diff))}")

    X, Y = np.meshgrid(x_i, t_j)
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, res1, cmap='inferno')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u(x,t)')
    plt.show()

    X, Y = np.meshgrid(x_i, t_j)
    fig = plt.figure(figsize=(6, 6), dpi=80)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, res2, cmap='inferno')
    ax.set_xlabel('x')
    ax.set_ylabel('t')
    ax.set_zlabel('u(x,t)')
    plt.show()
