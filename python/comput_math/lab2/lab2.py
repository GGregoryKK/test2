import math

from prettytable import PrettyTable
import matplotlib.pyplot as plt


def func(x: float, d: int = 1):
    return 1 / (1 + d * x * x)


def equally_spaced_nodes(a, b, n=4):
    if abs(a) == b:
        h = b / (n - 1)
    else:
        h = (b - a) / (n - 1)
    for i in range(2 * n - 1):
        yield a + i * h


def chebyshev_nodes(a, b, n):
    for i in range(n + 1):
        yield (a + b) / 2 + (b - a) / 2 * math.cos((2 * i - 1) * math.pi / (2 * n))


def lagr0(f, t, a, b, n, nodes):
    p = 0
    for j in nodes(a, b, n):
        p1 = 1
        p2 = 1
        for i in nodes(a, b, n):
            if i == j:
                p1 *= 1
                p2 *= 1
            else:
                p1 = p1 * (t - i)
                p2 = p2 * (j - i)
        p += f(t) * p1 / p2
    return p


def task_1():
    print(f"--" * 15, f"task 1.1", f"--" * 15)
    a, b, n = -5, 5, 4
    xx = [i for i in equally_spaced_nodes(a, b, n * 2 - 1)]
    interpolated = [lagr0(func, x, a, b, n, equally_spaced_nodes) for x in xx]
    table = PrettyTable()
    table.add_column("x", xx)
    table.add_column("y", [func(i) for i in xx])
    table.add_column("Li", interpolated)
    print(table)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")

    fig, ax = plt.subplots()
    plt.title("Lagrange (equally_spaced_nodes)")
    plt.plot(xx, interpolated, 'orange', label='interpolated')
    plt.scatter(xx, [func(i) for i in xx], label='$y_i$')
    xxx = [i for i in equally_spaced_nodes(a, b, n * 1000)]
    plt.plot(xxx, [func(i) for i in xxx], 'green', label='func')
    plt.legend(loc='best')
    plt.grid()
    # plt.show()
    plt.savefig("Lagrange_(equally_spaced_nodes).png")

    print(f"--" * 15, f"task 1.2", f"--" * 15)
    a, b, n = -5, 5, 13
    xx = [i for i in chebyshev_nodes(a, b, n)]
    xx.sort()
    interpolated = [lagr0(func, x, a, b, n, chebyshev_nodes) for x in xx]
    table = PrettyTable()
    table.add_column("x", xx)
    table.add_column("y", [func(i) for i in xx])
    table.add_column("Li", interpolated)
    print(table)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")

    fig, ax = plt.subplots()
    plt.title("Lagrange (chebyshev_nodes)")

    plt.plot(xx, interpolated, 'orange', label='interpolated')
    plt.scatter(xx, [func(i) for i in xx], label='$y_i$')
    xxx = [i for i in chebyshev_nodes(a, b, n * 1000)]
    plt.plot(xxx, [func(i) for i in xxx], 'green', label='func')

    plt.legend(loc='best')
    plt.grid()
    # plt.show()
    plt.savefig("Lagrange_(chebyshev_nodes).png")


def func_1(x):
    return abs(x)


def func_2(x):
    return x * x * x


def func_3(x):
    return math.cos(math.pi * x / 2)


def task_2(f, s1, s2):
    print(f"--" * 15, f"task 2.1 {s1}", f"--" * 15)
    a, b, n = -1, 1, 4
    xx = [i for i in equally_spaced_nodes(a, b, n * 2 - 1)]
    interpolated = [lagr0(f, x, a, b, n, equally_spaced_nodes) for x in xx]
    table = PrettyTable()
    table.add_column("x", xx)
    table.add_column("y", [f(i) for i in xx])
    table.add_column("Li", interpolated)
    print(table)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")

    fig, ax = plt.subplots()
    plt.title(f"Lagrange (equally_spaced_nodes) {s1}")
    plt.plot(xx, interpolated, 'orange', label='interpolated')
    plt.scatter(xx, [f(i) for i in xx], label='$y_i$')
    xxx = [i for i in equally_spaced_nodes(a, b, n * 1000)]
    plt.plot(xxx, [f(i) for i in xxx], 'green', label=f'{s2}')
    plt.legend(loc='best')
    plt.grid()
    # plt.show()
    plt.savefig(f"Lagrange_(equally_spaced_nodes)_{s1}.png")


if __name__ == "__main__":
    task_1()
    # task_2(func_1, "A", "f(x)=|x|")
    # task_2(func_2, "B", "$x^3$")
    # task_2(func_3, "C", "$cos(\\frac{\pi x}{2})$")
