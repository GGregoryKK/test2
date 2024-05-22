import copy


def dot_product(u, v):
    result = 0
    for i in range(len(u)):
        result += u[i] * v[i]
    return result


def simple_iteration(A, b, x0, max_iterations=1000, eps=1e-6):
    n = len(b)
    x = x0.copy()
    for i in range(max_iterations):
        x_new = x.copy()
        for j in range(n):
            s = sum(A[j][k] * x_new[k] for k in range(n) if k != j)  # (1.5)
            x_new[j] = (b[j] - s) / A[j][j]  # (1.5)
        if max(abs(x_new[i] - x[i]) for i in range(n)) < eps:  # (1.6)
            return x_new
        x = x_new
    return x


def speedy_descent(A, b, x0, max_iter=1000, tol=1e-8):
    r = [bi - dot_product(A[i], x0) for i, bi in enumerate(b)]  # (2.2)
    x = x0.copy()

    for i in range(max_iter):
        tau = dot_product(r, r) / dot_product(r, [dot_product(A[i], r) for i in range(len(A))])  # (2.5)
        x = [xi + tau * ri for xi, ri in zip(x, r)]  # x^{k+1} = tau^{k+1} r^{k}
        r_new = [bi - dot_product(A[i], x) for i, bi in enumerate(b)]  # (2.3) r^{k+1}
        if (dot_product(r_new, r_new)) ** .5 < tol:
            break
        r = r_new.copy()
    return x


def min_residuals(A, b, x0, max_iter=1000, eps=1e-8):
    r = [bi - dot_product(A[i], x0) for i, bi in enumerate(b)]  # (2.2)
    x = x0.copy()
    for i in range(max_iter):
        Ar = [dot_product(A[i], r) for i in range(len(A))]
        tau = dot_product(r, Ar) / dot_product(Ar, Ar)  # (2.4)
        x = [xi + tau * ri for xi, ri in zip(x, r)]  # x^{k+1} = tau^{k+1} r^{k}
        r_new = [bi - dot_product(A[i], x) for i, bi in enumerate(b)]  # (2.3) r^{k+1}
        if (dot_product(r_new, r_new)) ** .5 < eps:
            break
        r = r_new.copy()
    return x


from labs.lab7.lab7 import inverse

if __name__ == "__main__":
    A = [[4, 1, 2], [3, 5, 1], [1, 1, 6]]
    b = [4, 7, 3]
    x0 = [0, 0, 0]
    eps = 1e-5
    max_iter = 1000

    print(f"simple_iteration :\n{simple_iteration(copy.deepcopy(A), b, x0, max_iter, eps)}\n")

    print(f"speedy_descent :\n{speedy_descent(copy.deepcopy(A), b, x0, max_iter, eps)}\n")

    print(f"minimal_residual :\n{min_residuals(copy.deepcopy(A), b, x0, max_iter, eps)}\n")

    print(f"gaus:\n{inverse(A.copy(), b)}")
