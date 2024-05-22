import random
import copy
import math


def dot(x, y):
    return sum([xi * yi for xi, yi in zip(x, y)])


def norm(x):
    return sum([xi ** 2 for xi in x]) ** 0.5


def find_largest_eigenvalue(A, eps=1e-10):
    n = len(A)
    x = [random.random() for _ in range(n)]
    x_norm = norm(x)
    x = [xi / x_norm for xi in x]
    while True:
        y = [sum([A[i][j] * x[j] for j in range(n)]) for i in range(n)]  # (3)   x^{k+1}
        lmbda = dot(y, x) / dot(x, x)  # (2) lamda^{k}
        y_norm = norm(y)  # || x^{k+1} ||
        err = norm([yi - lmbda * xi for xi, yi in zip(x, y)]) / y_norm  # (4)

        if err < eps:
            break
        x = [yi / y_norm for yi in y]

    x_norm = norm(x)
    x = [xi / x_norm for xi in x]
    return lmbda, x


def find_smallest_eigenvalue(A, eps=1e-10):
    lmbda1, x1 = find_largest_eigenvalue(A, eps)
    n = len(A)
    B = [[A[i][j] - lmbda1 * (i == j) for j in range(n)] for i in range(n)]
    x = x1
    x_norm = norm(x)
    x = [xi / x_norm for xi in x]
    while True:
        y = [sum([B[i][j] * x[j] for j in range(n)]) for i in range(n)]     # (3)   x^{k+1}
        mu = dot(y, x) / dot(x, x)  # (2) lamda^{k}
        y_norm = norm(y)
        err = norm([yi - mu * xi for xi, yi in zip(x, y)]) / y_norm     # (4)
        if err < eps:
            break
        x = [yi / y_norm for yi in y]
    x_norm = norm(x)
    x = [xi / x_norm for xi in x]
    return lmbda1 + mu, x


def jacobi_rotation(A, eps=1e-10):
    n = len(A)
    eig_vecs = [[0.0] * n for i in range(n)]
    for i in range(n):
        eig_vecs[i][i] = 1.0

    while True:
        # Находим наибольший недиагональный элемент матрицы.
        p, q = 0, 0
        max_off_diag = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > max_off_diag:
                    max_off_diag = abs(A[i][j])
                    p, q = i, j
                    print(f"k = {p}; l = {q}")

        if max_off_diag < eps:
            # Матрица уже диагональна.
            eig_vals = [A[i][i] for i in range(n)]
            return eig_vals, eig_vecs

        # Вычисляем угол вращения.
        app = A[p][p]
        aqq = A[q][q]
        apq = A[p][q]
        phi = 0.5 * math.atan2(2.0 * apq, aqq - app)

        # Вычисляем элементы матрицы вращения.
        c = math.cos(phi)
        s = math.sin(phi)

        # Выполняем вращение.
        for i in range(n):
            if i != p and i != q:
                api = A[p][i]
                aqi = A[q][i]
                A[p][i] = c * api - s * aqi
                A[q][i] = s * api + c * aqi
                A[i][p] = A[p][i]
                A[i][q] = A[q][i]

            # Обновляем собственные векторы.
            epi = eig_vecs[i][p]
            eqi = eig_vecs[i][q]
            eig_vecs[i][p] = c * epi - s * eqi
            eig_vecs[i][q] = s * epi + c * eqi

        # Обновляем диагональные элементы.
        app = A[p][p]
        aqq = A[q][q]
        apq = A[p][q]
        A[p][p] = c * c * app - 2.0 * c * s * apq + s * s * aqq
        A[q][q] = s * s * app + 2.0 * c * s * apq + c * c * aqq
        A[p][q] = 0.0
        A[q][p] = 0.0


if __name__ == "__main__":
    A = [[3, 9],
         [9, 4]]

    print(f"max lambda {find_largest_eigenvalue(copy.deepcopy(A))}")
    print(f"min lambda {find_smallest_eigenvalue(copy.deepcopy(A))}")

    l, al = jacobi_rotation(copy.deepcopy(A))
    print(f"\n\n{l}\n")
    for i in al:
        print(f"vec {i}\n")
