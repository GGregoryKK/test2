from typing import Callable

from matplotlib import pyplot as plt
import numpy as np

N = 100
k0 = 66


def func(x: float) -> float:
    return x * (np.pi - x) * (np.pi - x)


def noize(x: float, n: int = N) -> float:
    t = (x - (k0 + 0.5) * np.pi / n)
    return 1 / (np.log(n)) * np.cos(n * x) * t / abs(t)


def noize_f(x: float, n: int = N) -> float:
    return func(x) + noize(x=x, n=n)


def sinc(f: Callable[[float], float], x: float, n: int = N):
    return sum(
        [0 if x == 0 else ((-1) ** k * np.sin(n * x) / (n * x - k * np.pi)) * f(k * np.pi / n) for k in range(n)])


def my_operators(f: Callable[[float], float], x: float, n: int = N):
    summ = 0
    for k in range(n):
        d0, d1 = (N * x - k * np.pi) * (-1) ** (k), (N * x - (k + 1) * np.pi) * (-1) ** (k + 1)
        if d0 == 0 or d1 == 0:
            continue
        d = np.sin(N * x)
        s0 = d / d0
        s1 = d / d1
        summ += (s0 + s1) * f(k * np.pi / N)
    return summ / 2


if __name__ == "__main__":
    x = np.linspace(0, np.pi, N)

    res_default_sinc = [sinc(noize_f, xi, N) for xi in x]
    res_op_sinc = [my_operators(noize_f, xi, N) for xi in x]

    plt.plot(x, [func(i) for i in x], label="Func", linewidth=2)
    plt.plot(x, [noize_f(xi, N) for xi in x], label="Func with noize", linewidth=2)
    plt.plot(x, res_default_sinc, label="Sinс", linewidth=1.5)
    plt.plot(x, res_op_sinc, label="$B_{\\lambda}$", linewidth=2)
    plt.grid()
    plt.legend()
    plt.show()
