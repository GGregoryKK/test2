# https://ru.wikipedia.org/wiki/%D0%9E%D1%81%D1%86%D0%B8%D0%BB%D0%BB%D1%8F%D1%82%D0%BE%D1%80_%D0%92%D0%B0%D0%BD_%D0%B4%D0%B5%D1%80_%D0%9F%D0%BE%D0%BB%D1%8F
import numpy as np
import matplotlib.pyplot as plt


def plot1(t, u):
    plt.figure()
    x, y = [i[0] for i in u], [i[1] for i in u]
    plt.subplot(2, 2, 1)
    plt.plot(t, x)
    plt.xlabel('t')
    plt.ylabel('x')
    plt.grid()
    plt.title('coordinate')

    plt.subplot(2, 2, 2)
    plt.plot(t, y)
    plt.xlabel('t')
    plt.ylabel('y')
    plt.grid()
    plt.title('velocity')

    plt.subplot(2, 2, 3)
    plt.plot(x, y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    plt.title('phase picture')

    plt.subplot(2, 2, 4)
    plt.plot(t, u)
    plt.xlabel('t')
    plt.ylabel('y and x')
    plt.title('velocity')
    plt.grid()
    plt.tight_layout()

    plt.show()


def f(y, t):
    mu = .5
    f2 = y[0] / mu
    f1 = mu * (f2 - f2 ** 3 / 3 - y[1])

    return np.array([f1, f2])


def ff(y, t):
    mu = .5
    f1 = mu * (y[0] - y[0] ** 3 / 3 - y[1])
    f2 = y[0] / mu
    return [f1, f2]


def odeint(f, X_0, steps):
    usol = [X_0]
    dt = np.diff(steps)[0]
    u = np.copy(X_0)
    for t in steps[1:]:
        u1 = f(u + 0.5 * dt * f(u, t), t + 0.5 * dt)
        u2 = f(u + 0.5 * dt * u1, t + 0.5 * dt)
        u3 = f(u + dt * u2, t + dt)
        u = u + (1 / 6) * dt * (f(u, t) + 2 * u1 + 2 * u2 + u3)
        usol.append(u)
    return usol


def odeint_euler_mod(f, X_0, steps):
    usol = [X_0]
    dt = np.diff(steps)[0]
    u = np.copy(X_0)
    for t in steps[1:]:
        u1 = u + dt * f(u, t + dt)
        u = u + .5 * dt * (f(u, t + dt) + f(u1, t + dt))
        usol.append(u)
    return usol


def wut(f, X_0, steps, eps):
    p_len = len(steps)
    solve = odeint_euler_mod(f, X_0, steps)
    p_solve_l = [i[0] for i in solve]
    p_solve_r = [i[1] for i in solve]
    max_div_l = 0
    max_div_r = 0
    left, right = False, False
    del solve
    while True:
        new_len = p_len * 2
        t = np.linspace(steps[0], steps[len(steps) - 1], new_len)
        new_solve = odeint(f, X_0, t)
        if not left:
            new_solve_l = [i[0] for i in new_solve]
            for i in range(1, int(new_len / 2)):
                m = abs(p_solve_l[i] - new_solve_l[i * 2])
                if m > max_div_l:
                    max_div_l = m
            if max_div_l < eps:
                left = True
        if not right:
            new_solve_r = [i[1] for i in new_solve]
            for i in range(1, int(new_len / 2)):
                m = abs(p_solve_r[i] - new_solve_r[i * 2])
                if m > max_div_r:
                    max_div_r = m
            if max_div_r < eps:
                right = True

        if right and left:
            return new_solve, t
        else:
            max_div_r, max_div_l = 0, 0
            p_len = new_len
            p_solve_l = [i[0] for i in new_solve]
            p_solve_r = [i[1] for i in new_solve]


if __name__ == "__main__":
    to = 0
    t_end = 20

    tt = np.linspace(to, t_end, 3)
    start_yo = np.array([0, 1])

    u, t = wut(f, start_yo, tt, .001)

    plot1(t, u)
