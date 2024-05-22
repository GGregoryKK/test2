import numpy as np


def phi(theta: float) -> float:
    return np.sin(theta * np.pi)


def iphi(ksi: float, omega1: int | float, omega2: int | float) -> float:
    a = omega2 / (omega2 - omega1)
    if 0 <= ksi < a:
        return omega2 * phi(ksi / a)    # (7)
    elif a <= ksi <= 1:
        return omega1 * phi((1 - ksi) / (1 - a))    # (7)
    else:
        return 0    # мое условие)


def farg(alpha: float) -> float:
    return np.sign(alpha) * abs(int(alpha) - alpha) # дробная часть


def u_xt(u_lst: np.ndarray, x_lst: np.ndarray, t_lst: np.ndarray, omega: tuple):
    xi = 0
    d_omega = omega[1] - omega[0]
    for x in x_lst:
        ti = 0
        for t in t_lst:
            arg1 = farg((t + omega[1] * x) / d_omega)
            arg2 = farg((t + omega[0] * x) / d_omega)
            u_lst[ti][xi] = (iphi(arg1, omega[0], omega[1]) - iphi(arg2, omega[0], omega[1])) / d_omega     # (6)
            ti += 1
        xi += 1
    return u_lst
