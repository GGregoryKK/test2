from typing import Callable


def sign(a: float):
    if a > 0:
        return 1
    elif a < 0:
        return -1
    else:
        return 0


def ping(a: float, b: float, f: Callable[[float], float]):
    if sign(f(a)) == sign(f(b)):
        return False
    else:
        return True


def half_divide_method(a: float | int, b: float | int, f: Callable[[float], float], eps: float):
    h = eps * b * 2 if abs(a) == b else eps * (b - a)
    x = a
    l = []
    while x < b:
        k = x
        x += h
        if ping(x, k, f):
            left = k
            right = x
            p = 0
            while abs(right - left) >= eps or (f(right) >= eps) or abs(right - left) / right >= eps:
                p = (left + right) / 2
                if ping(left, p, f):
                    right = p
                else:
                    left = p
            l.append(p)
    return l


def newton(f: Callable[[float], float], f_prime: Callable[[float], float], x0: float,
           eps: float = 1e-7, kmax: int = 1e3) -> float:
    x, x_prev, i = x0, x0 + 2 * eps, 0

    while abs(x - x_prev) >= eps and i < kmax:
        x, x_prev, i = x - f(x) / f_prime(x), x, i + 1

    return x


def func(x: int | float):
    """
    x ^ 2 - 4
    """
    return x ** 2 - 4


def dx_func(x: int | float):
    return 2 * x


if __name__ == "__main__":
    a, b = -5, 5
    n = 100
    eps = 1e-3

    h_res = half_divide_method(a, b, func, eps)
    print(f"Root of the function: {func.__doc__} at [-5, 5].\n")
    for i in h_res:
        print(f"HDM {round(i, 5)}")
        print(f"Newton {round(newton(func, dx_func, i, eps=eps), 5)}\n")
