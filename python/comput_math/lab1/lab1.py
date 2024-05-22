from math import exp, factorial, sin, cos, pow, radians
from prettytable import PrettyTable


def _task_1_1(ek: int | float) -> (float, int):
    count = 0
    while True:
        if ek / 2 == 0:
            return ek, count
        count += 1
        ek /= 2


def _task_1_2(ek: int | float) -> (float, int):
    count = 0
    while 1 + ek > 1:
        count += 1
        ek /= 2
    return ek, count


def _task_2_1(x: int) -> (float, int):
    sum = 1
    fact = 1
    result = 1
    c = 1
    while 1 + sum > 1:
        fact *= c
        sum = (x ** c / fact)
        result += sum
        c += 1
    return result


def _task_2_3_cos(x: int | float) -> (float, int):
    sum = 1
    result = 0
    c = 0
    while 1 + abs(sum) > 1:
        sum = (pow(x, 2 * c) / factorial(2 * c)) * (-1) ** (c & 1)
        result += sum
        c += 1
    return result


def _task_2_3_sin(x: int | float) -> (float, int):
    sum = 1
    result = 0
    c = 0
    while 1 + abs(sum) > 1:
        sum = (pow(x, 2 * c + 1) / factorial(2 * c + 1)) * (-1) ** (c & 1)
        result += sum
        c += 1
    return result


def _task_3_1(n: int, eps: float) -> float:
    e = exp(-1) + eps
    for k in range(2, n + 1):
        e = 1 - k * e
    return e


def _task_3_2(n: int, eps: float) -> float:
    e = eps
    for k in range(n, 0, -1):
        e = (1 - e) / k
    return e


def task_1_1() -> None:
    print(f"--" * 15, f"task 1.1", f"--" * 15)
    number = 1
    ek, k = _task_1_1(ek=number)
    print(f"ek in N; num = {number}:\n\tk = {k};\tek = {ek}")

    number = 1.35
    ek, k = _task_1_1(ek=number)
    print(f"ek in R; num = {number}:\n\tk = {k};\tek = {ek}")
    print(f"--" * 15 + f"-" * 11 + f"--" * 15)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_1_2() -> None:
    print(f"--" * 15, f"task 1.2", f"--" * 15)

    number = 1
    fk, k = _task_1_2(ek=number)
    print(f"fk in N; num = {number}:\n\tk = {k};\tfk = {fk}")

    number = 1.35
    fk, k = _task_1_2(ek=number)
    print(f"fk in R; num = {number}:\n\tk = {k};\tfk = {fk}")
    print(f"--" * 15 + f"-" * 11 + f"--" * 15)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_2_1() -> None:
    print(f"--" * 15, f"task 2.1", f"--" * 15)
    x = [(_task_2_1(10 * i), exp(10 * i), i * 10) for i in range(5)]
    for res in x:
        print(f"x = {res[2]}:\n\tresult: {res[0]};\n\tstd: {res[1]};\n\tdelta {abs(res[0] - res[1])}.\n")
    print(f"--" * 15 + f"-" * 11 + f"--" * 15)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_2_2() -> None:
    print(f"--" * 15, f"task 2.2", f"--" * 15)
    x = [(1 / _task_2_1(10 * i), exp(-10 * i), i * 10) for i in range(5)]
    for res in x:
        print(f"x = {-res[2]}:\n\tresult: {res[0]};\n\tstd: {res[1]};\n\tdelta {abs(res[0] - res[1])}.\n")
    print(f"--" * 15 + f"-" * 11 + f"--" * 15)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_2_3() -> None:
    print(f"--" * 15, f"task 2.3", f"--" * 15)
    x_cos = [(_task_2_3_cos(10 * i), cos(10 * i), i * 10) for i in range(5)]
    x_sin = [(_task_2_3_sin(10 * i), sin(10 * i), i * 10) for i in range(5)]
    for r_cos, r_sin in zip(x_cos, x_sin):
        print(
            f"x = {r_cos[2]}:\n\tresult cos: {r_cos[0]};\tresult sin: {r_sin[0]};\n\tstd cos: {r_cos[1]};\tstd sin: {r_sin[1]};\n\tdelta cos: {abs(r_cos[0] - r_cos[1])};\tdelta sin: {abs(r_sin[0] - r_sin[1])}.\n")
    print(f"--" * 15 + f"Radians" + f"--" * 15)
    x_cos = [(_task_2_3_cos(radians(10 * i)), cos(radians(10 * i)), radians(i * 10)) for i in range(5)]
    x_sin = [(_task_2_3_sin(radians(10 * i)), sin(radians(10 * i)), radians(10 * i)) for i in range(5)]
    for r_cos, r_sin in zip(x_cos, x_sin):
        print(
            f"x = {r_cos[2]}:\n\tresult cos: {r_cos[0]};\tresult sin: {r_sin[0]};\n\tstd cos: {r_cos[1]};\tstd sin: {r_sin[1]};\n\tdelta cos: {abs(r_cos[0] - r_cos[1])};\tdelta sin: {abs(r_sin[0] - r_sin[1])}.\n")
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_3_1() -> None:
    print(f"--" * 15, f"task 3.1", f"--" * 15)
    epsilon = [0, 1e-7, 1e-6, 1e-5]
    fct20 = factorial(20)
    for eps in epsilon:
        ek = _task_3_1(20, eps)
        err = fct20 * eps
        print(f" k = 20; ek = {ek}; k! * eps = {err}; ek + k! * eps = {ek + err}")
        print(f"--" * 15 + f"-" * 11 + f"--" * 15)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_3_2() -> None:
    print(f"--" * 15, f"task 3.2", f"--" * 15)
    epsilon = [0, 1e-2, 1e-1, 1]
    fct20 = factorial(20)
    for eps in epsilon:
        ek = _task_3_2(30, eps)
        err = fct20 * eps
        print(f" k = 20; ek = {ek}; k! * eps = {err}; ek + k! * eps = {ek + err}")
        print(f"--" * 15 + f"-" * 11 + f"--" * 15)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


def task_4() -> None:
    def _li(i):
        prod = 1

        for j in range(1, 21):
            if i == j:
                pass
            else:
                prod *= (i - j)

        return i ** 19 / prod

    table = PrettyTable()
    roots = [i for i in range(1, 21)]
    li = [_li(root) for root in roots]
    a = 1e-7
    delta_x = [l * a for l in li]
    i_plus_delta_x = [i + x for i, x in zip(roots, delta_x)]
    # table.field_names = ["i", "Li", "delta x", "i + delta x"]
    table.add_column("i", roots)
    table.add_column("Li", li)
    table.add_column("delta x", delta_x)
    table.add_column("i + delta x", i_plus_delta_x)
    print(f"--" * 15, f"task 4", f"--" * 15)
    print(table)
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "")
    print(f"--" * 15 + f"-" * 11 + f"--" * 15, "\n\n")


if __name__ == "__main__":
    task_1_1()
    task_1_2()
    task_2_1()
    task_2_2()
    task_2_3()
    task_3_1()
    task_3_2()
    task_4()
    pass
