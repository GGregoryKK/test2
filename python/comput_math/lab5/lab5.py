import random

import matplotlib.pyplot as plt
from numpy import log, exp, array, linalg


def lsm(x: list, y: list, deg: int = 2, power: bool = False):
    sX = 0
    sY = 0
    sXY = 0
    sXX = 0
    n = len(y)
    if power:
        for i in range(len(x)):
            sX += log(x[i])
            sXX += log(x[i]) ** 2
            if y[i] > 0:
                sY += log(y[i])
                sXY += log(x[i]) * log(y[i])

        b = (n * sXY - sX * sY) / (n * sXX - sX ** 2)
        # ln_a = (sY - sX * b) / n
        a = exp((sY - b * sX) / n)
        return [a * i ** b for i in x]
    else:
        powered_x_sums = array([sum(x ** i) for i in range(len(x))])
        a = array([array([0 for _ in range(deg)]) for _ in range(deg)])
        for i in range(deg):
            for j in range(deg):
                a[i][j] = powered_x_sums[i + j]
        inv_a = linalg.inv(a)
        b = [sum(y), sum(x * y), sum(x ** 2 * y)]
        b = [b[i] for i in range(deg)]
        xx = inv_a.dot(b)
        if deg == 2:
            return [xx[0] + xx[1] * i for i in x]
        if deg == 3:
            return [xx[0] + xx[1] * i + xx[2] * i ** 2 for i in x]


def plot(x: list, y: list, y1: list, title: str ="", show: bool = True, save: bool = False) -> None:
    fig, ax = plt.subplots()
    plt.title(f"{title}")
    plt.scatter(x, y, label='source data')
    plt.scatter(x, y1, label='$LSM$')
    plt.legend(loc='best')
    plt.grid()
    if show:
        plt.show()
    if save:
        title = title.replace(' ', '_').replace('\\alpha', 'a').replace('$', '').replace("+", "").replace("^", "")
        plt.savefig(f"{title}.png")


if __name__ == "__main__":
    n = 100
    x = [random.randint(1, 200) for _ in range(n)]
    y = [random.randint(1, 400) for _ in range(n)]
    z = lsm(array(x), array(y), power=True)
    plot(x, y, z, title="$\\alpha x^b$", show=False, save=True)
    z = lsm(array(x), array(y), 2)
    plot(x, y, z, title="$\\alpha_0 + \\alpha_1 x$", show=False, save=True)
    z = lsm(array(x), array(y), 3)
    plot(x, y, z, title="$\\alpha_0 + \\alpha_1 x + \\alpha_2 x^2$", show=False, save=True)
