import numpy as np
from math import sqrt

from lab2.lab2 import equally_spaced_nodes


def plot(x: list, y: list, y1: list, title: str = "", show: bool = True, save: bool = False) -> None:
    fig, ax = plt.subplots()
    plt.title(f"{title}")
    plt.plot(x, y, label='func')
    plt.plot(x, y1, label='S')

    plt.legend(loc='best')
    plt.grid()
    if show:
        plt.show()
    if save:
        title = title.replace(' ', "").replace("=", "_")
        plt.savefig(f"{title}.png")


def cubic_interp1d(x0, x, y):
    x = np.asfarray(x)
    y = np.asfarray(y)

    if np.any(np.diff(x) < 0):
        indexes = np.argsort(x)
        x = x[indexes]
        y = y[indexes]

    size = len(x)

    xdiff = np.diff(x)
    ydiff = np.diff(y)

    li = np.empty(size)
    li_1 = np.empty(size - 1)
    z = np.empty(size)

    li[0] = sqrt(2 * xdiff[0])
    li_1[0] = 0.0
    B0 = 0.0
    z[0] = B0 / li[0]

    for i in range(1, size - 1, 1):
        li_1[i] = xdiff[i - 1] / li[i - 1]
        li[i] = sqrt(2 * (xdiff[i - 1] + xdiff[i]) - li_1[i - 1] * li_1[i - 1])
        bi = 6 * (ydiff[i] / xdiff[i] - ydiff[i - 1] / xdiff[i - 1])
        z[i] = (bi - li_1[i - 1] * z[i - 1]) / li[i]

    i = size - 1
    li_1[i - 1] = xdiff[-1] / li[i - 1]
    li[i] = sqrt(2 * xdiff[-1] - li_1[i - 1] * li_1[i - 1])
    bi = 0.0
    z[i] = (bi - li_1[i - 1] * z[i - 1]) / li[i]

    i = size - 1
    z[i] = z[i] / li[i]
    for i in range(size - 2, -1, -1):
        z[i] = (z[i] - li_1[i - 1] * z[i + 1]) / li[i]

    index = x.searchsorted(x0)
    np.clip(index, 1, size - 1, index)

    xi1, xi0 = x[index], x[index - 1]
    yi1, yi0 = y[index], y[index - 1]
    zi1, zi0 = z[index], z[index - 1]
    hi1 = xi1 - xi0

    f0 = zi0 / (6 * hi1) * (xi1 - x0) ** 3 + \
         zi1 / (6 * hi1) * (x0 - xi0) ** 3 + \
         (yi1 / hi1 - zi1 * hi1 / 6) * (x0 - xi0) + \
         (yi0 / hi1 - zi0 * hi1 / 6) * (xi1 - x0)
    return f0


def func(x, d=21):
    return 1 / (1 + d * x * x)


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    n_points = 6
    a, b = -1, 1
    for n_points in range(2, 6):
        x = [i for i in equally_spaced_nodes(-1, 1, n_points)]
        y = [func(i) for i in x]
        x_new = [i for i in equally_spaced_nodes(-1, 1, n_points * 3)]
        yy = [func(i) for i in x_new]
        plot(x_new, yy, cubic_interp1d(x_new, x, y), title=f"n_points = {n_points}", show=False, save=True)
