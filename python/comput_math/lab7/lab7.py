def appI(A):  # Расширяет матрицу А (Приписывает единичную справа.) A(nn) -> A(n2n)
    l = len(A)
    i = [0 for i in range(l)]
    for j in range(l):
        i[j] = 1
        A[j].extend(i)
        i[j] = 0
    return A


def operators(ch, a, b):  # позволяет делать линейные преобразования (Для метода Гаусса нах обр матрицы)
    c = []
    if ch == "-":
        for i in range(len(a)):
            c.append(a[i] - b[i])
        return c

    elif ch == "+":
        for i in range(len(a)):
            c.append(a[i] + b[i])
        return c

    elif ch == "*":
        for i in a:
            c.append(i * b)
        return c

    elif ch == "/":
        for i in a:
            c.append(i / b)
        return c


def E(n, bool):  # Создает единичную (bool == 1) или нулевую (bool == 0) матрицу размерности nxn
    i = [0 for i in range(n)]
    E = []
    if bool == 1:
        for j in range(n):
            i[j] = 1
            E.append(list(i))
            i[j] = 0
    else:
        for j in range(n):
            E.append(list(i))
    return E


def inverse(M, inverse=True):  # Метод Гаусса
    if isinstance(inverse, bool):
        n = len(M)
        M = appI(M)
    else:
        n = 1
        for i in range(len(inverse)):
            M[i].append(inverse[i])

    l = len(M)
    for i in range(l):
        z = 0
        if M[i][i] == 0:  # Если элемент на гравной диагонали == 0, то меняем эту строку на любую другую, где ai=j != 0
            for h in range(l):
                if M[h][i] != 0:
                    z = h
        if z != 0:  # Прямой ход
            M[i], M[z] = M[z], M[i]
            M[i] = operators('/', M[i], M[i][i])
            j = i + 1
            for k in range(j, l):
                c = M[k][i] / M[i][i]
                M[k] = operators('-', M[k], operators('*', M[i], c))
        else:
            M[i] = operators('/', M[i], M[i][i])
            j = i + 1
            for k in range(j, l):
                c = M[k][i] / M[i][i]
                M[k] = operators('-', M[k], operators('*', M[i], c))
    l1 = l - 1
    for i in range(l1, 0, -1):  # Обратный ход
        M[i] = operators('/', M[i], M[i][i])
        j = i - 1
        for k in range(j, -1, -1):
            c = M[k][i] / M[i][i]
            M[k] = operators('-', M[k], operators('*', M[i], c))

    if isinstance(inverse, bool):   # Отделение обратной от расширенной матрицы
        N = E(len(M), 0)
        for i in range(len(M)):
            for j in range(n):
                N[i][j] = M[i][len(M) + j]
    else:
        N = []
        for i in range(len(M)):
            for j in range(n):
                N.append(M[i][len(M) + j])
    return N


if __name__ == "__main__":
    a = [
            [1, 2],
            [1, 3]
        ]

    b = [1, 2]
    print(f"SLAU:\n{inverse(a.copy(), b)}")

    a = [
            [7, 4],
            [2, 1]
        ]
    print(f"Inverse:\n{inverse(a)}")
