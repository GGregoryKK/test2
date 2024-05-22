from math import sin, cos, cosh, sinh


def func1(x, a=2):
    return sin(a * x) / x


def func2(x):
    return (cos(x) - 1) / x


def func3(x, a=2):
    return sinh(a * x) / x


def func4(x, a=2):
    return (1 - cos(a * x)) / x ** 2


def func5(x):
    return (cosh(x) - 1) / x


def func6(x, a=2):
    return (1 - cosh(a * x)) / x ** 2


def simp(f, a, b, n=100):
    s = f(a) + f(b)
    h = (b - a) / (n * 2)
    for i in range(1, n):
        s += 4 * f(2 * i * h - 1 + a) + 2 * f(2 * i * h + a)
    return h * (s + 4 * f(2 * n * h - 1 + a)) / 3


def trap(f, a, b, n=100):
    s = 0
    h = (b - a) / n
    s = (f(a) + f(b)) / 2
    for i in range(1, n):
        s += f(a + i * h)
    return h * s


def req(f, a, b, n=100):
    s = 0
    h = (b - a) / n
    for i in range(n):
        s += f((a + i * h + h / 2))
    return h * s


def hi(n: int):
    for i in range(n):
        pass


if __name__ == "__main__":
    a, b = 1e-100, 1
    methds = {"simp": simp, "req": req, "trap": trap}
    funcs = {"sin(a * x) / x": func1, "(cos(x) - 1) / x": func2, "sinh(a * x) / x": func3,
             "(1 - cos(a * x)) / x ** 2": func4, "(cosh(x) - 1) / x": func5, "(1 - cosh( a * x)) / x ** 2": func6}
    for name, func in funcs.items():
        print(f"{name}:")
        for k, v in methds.items():
            print(f"{k}: {v(func, a, b, 1000)}\n")
        print(f"\n\n")
