import matplotlib.pyplot as plt

from lab2.lab2 import equally_spaced_nodes


def plot(x: list, y: list, y1: list, y2: list, title: str = "", show: bool = True, save: bool = False) -> None:
    fig, ax = plt.subplots()
    plt.title(f"{title}")
    plt.plot(x, y, label='func')
    plt.plot(x, y1, label='Hermit y`(-1) = 1, y`(-1) = 1')
    plt.plot(x, y2, label='Hermit y`(-1) = 1, y`(-1) = -1')
    plt.legend(loc='best')
    plt.grid()
    if show:
        plt.show()
    if save:
        title = title.replace(' ', "").replace("=", "_")
        plt.savefig(f"{title}.png")


def func(x, d=5):
    return 1 / (1 + d * x * x)


def hermit(f, a, b, nu, d=21):
    h = 2
    x, y, z = [], [], []

    for i in equally_spaced_nodes(a, b, 11):
        x.append(i)
        y.append(f(i, d))
        k = (i - a) / h
        z.append((1 - 3 * k ** 2 + 2 * k ** 3) * nu[0] + (3 * k ** 2 - 2 * k ** 3) * nu[1] + \
                 nu[2] * h * (k - 2 * k ** 2 + k ** 3) - nu[3] * h * (k ** 2 - k ** 3))

    return x, y, z


if __name__ == "__main__":
    # task 2_2
    a, b = -1, 1
    for d in range(21):
        x, y, z = hermit(func, a, b, [func(-1, d), func(1, d), 1, 1])
        x, y, z1 = hermit(func, a, b, [func(-1, d), func(1, d), 1, -1])
        plot(x, y, z, z1, title=f"d = {d}", show=False, save=True)
