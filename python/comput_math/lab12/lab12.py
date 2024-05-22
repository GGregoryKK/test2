import numpy as np
from matplotlib import pyplot as plt
import math


def show_res(n_begin, n_end, h_begin, h_end, eps, m, x, y, yn, count):
    maxYn = 0
    print(f"ODU: \n {p.__doc__} y'' + {q.__doc__} y' + {r.__doc__} y = {f.__doc__}")
    print(f'Conditions: y({x[0]}) = {y[0]}, y({x[-1]}) = {y[-1]};')
    for i in range(n_end):
        d = abs(exact_solve(x[i]) - y[i])
        if d > maxYn:
            maxYn = d
    print(
        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print(
        '| n - begin | n - end | h - begin |     h - end    |  epsilon  | max differences Y | max differences Runge | count |')
    print(
        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print(
        '| {: ^9d} | {: ^7d} | {: ^9.3f} | {: ^14.6f} | {: ^9.5f} | {: ^17.5f} | {: ^21.5f} | {: ^5d} |'.format(n_begin,
                                                                                                                n_end,
                                                                                                                h_begin,
                                                                                                                h_end,
                                                                                                                eps,
                                                                                                                maxYn,
                                                                                                                m,
                                                                                                                count))
    print(
        '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Table:')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('|     X     |       Y       |       Y~      |     |Y - Y~|     |      Runge      |')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    for i in range(n_end):
        if i % 2 == 0:
            j = int(i / 2)
            print('| {: ^9.4f} | {: ^13.9f} | {: ^13.9f} | {: ^16.14f} | {: ^15.12f} |'.format(x[i], exact_solve(x[i]),
                                                                                               y[i],
                                                                                               abs(exact_solve(x[i]) -
                                                                                                   y[i]),
                                                                                               abs(yn[j] - y[i])))
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        else:
            print('| {: ^9.4f} | {: ^13.9f} | {: ^13.9f} | {: ^16.14f} | {: ^15s} |'.format(x[i], exact_solve(x[i]),
                                                                                            y[i],
                                                                                            abs(exact_solve(x[i]) -
                                                                                                y[i]),
                                                                                            'NO'))
            print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


def exact_solve(x):
    return (-2 * math.e * x - 2 * x + math.exp(1 - x) + math.exp(x)) / (
            1 + math.e)


def p(x):
    """ 1 """
    return 1


def q(x):
    """ + 0 """
    return 0


def r(x):
    """ -1 """
    return -1


def f(x):
    """ 2x """
    return 2 * x


def solve_pde(p, q, r, f, a, b, n, y0, yk):
    """
    Решение уравнения вида p(x) * y'' + q(x) * y' + r(x) * y = f(x) методом прогонки.
    p, q, r, f - функции, задающие коэффициенты уравнения.
    a, b - границы области определения уравнения.
    n - число узлов сетки.
    """
    h = (b - a) / (n - 1)
    x = [a + i * h for i in range(n)]

    a_vals = [0] * (n - 2)
    b_vals = [0] * (n - 2)
    c_vals = [0] * (n - 2)
    d_vals = [0] * (n - 2)

    for i in range(n - 2):
        a_vals[i] = p(x[i]) / h ** 2
        b_vals[i] = -(p(x[i]) + p(x[i])) / h ** 2 - q(x[i + 1]) / h
        c_vals[i] = p(x[i + 1]) / h ** 2 + q(x[i + 1]) / h
        d_vals[i] = f(x[i + 1]) - r(x[i + 1]) * h / 2

    d_vals[0] = d_vals[0] - a_vals[0] * y0
    d_vals[-1] = d_vals[-1] - c_vals[-1] * yk

    # Решение системы линейных уравнений методом прогонки
    y_vals = np.zeros(n)
    # Задание граничных условий
    y_vals[0] = y0
    y_vals[-1] = yk

    # Решение системы линейных уравнений методом прогонки
    y_vals[1:-1] = tridiagonal_matrix_algorithm(a_vals, b_vals, c_vals, d_vals)

    return x, y_vals


def tridiagonal_matrix_algorithm(a, b, c, d):
    """
    Реализация метода прогонки для решения трехдиагональных систем линейных уравнений.
    a, b, c - коэффициенты трехдиагональной матрицы системы.
    d - правая часть системы.
    """
    n = len(d)
    for i in range(1, n):
        m = a[i] / b[i - 1]
        b[i] = b[i] - m * c[i - 1]
        d[i] = d[i] - m * d[i - 1]

    x = np.zeros_like(d)
    x[-1] = d[-1] / b[-1]

    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]

    return x


def runge(p, q, r, f, a, b, n, y0, yk, eps):
    m = 0
    n2 = 2 * n
    c = 0
    while True:
        c += 1
        x1, y1 = solve_pde(p, q, r, f, a, b, n, y0, yk)
        x2, y2 = solve_pde(p, q, r, f, a, b, n2, y0, yk)
        for i in range(1, int(n / 2)):
            d = abs(y1[i] - y2[2 * i])
            if d > m:
                m = d
        if m < eps:
            break
        else:
            n = n2
            n2 *= 2
            m = 0

    return x2, y2, y1, c, n2, m


# ---------------------------------

def f1(t, y):
    dphi, dz = y
    return [dz, -9.81 / 10 * math.sin(dphi)]


def shoot(w):
    # начальные условия
    y0 = [0, w]
    # шаг
    h = .001
    # конечное время
    t_end = 1
    # количество шагов
    n = int(t_end / h)
    # метод Рунге-Кутты
    for i in range(n):
        k1 = f1(i * h, y0)
        k2 = f1(i * h + h / 2, [y0[j] + k1[j] / 2 for j in range(2)])
        k3 = f1(i * h + h / 2, [y0[j] + k2[j] / 2 for j in range(2)])
        k4 = f1(i * h + h, [y0[j] + k3[j] for j in range(2)])
        y0 = [y0[j] + h * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6 for j in range(2)]
    # возвращаем значение d phi / dt
    return y0[0]


def bisection(a, b, eps):
    """
    Функция, которая реализует метод бисекции для нахождения корня
    """
    # начальное приближение
    c = (a + b) / 2
    # значение функции в точке a
    fa = shoot(a)
    # значение функции в точке c
    fc = shoot(c)
    # пока не достигнуто требуемое значение точности
    while abs(fc) > eps:
        # если знак fc и fa совпадают
        if fc * fa > 0:
            # сдвигаем левую границу
            a = c
            fa = fc
        # иначе
        else:
            # сдвигаем правую границу
            b = c
        # вычисляем новое значение c
        c = (a + b) / 2
        fc = shoot(c)
    # возвращаем значение c
    return c


def task1():
    # начальное приближение для метода бисекции
    a = 0
    b = 2
    # требуемая точность
    eps = 1e-6
    # находим корень
    w = bisection(a, b, eps)
    # выводим результат
    print(f"Task1: \nw = {w}\n\n")


def task2():
    n = 4
    a, b = 0, 1
    x0, y0 = 0, 1
    xk, yk = 1, -1
    eps = 1e-1
    count = 0

    x, y, y1, c, n2, m = runge(p, q, r, f, a, b, n, y0, yk, eps)

    show_res(n, n2, (b - a) / n, (b - a) / n2, eps, m, x, y, y1, c)

    plt.plot(x, [exact_solve(xi) for xi in x], label="Exact")
    plt.plot(x, y, label="Runge run")
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    task1()
    task2()
