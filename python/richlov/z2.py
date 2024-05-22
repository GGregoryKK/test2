import numpy as np


def phi(theta: float) -> float:
    if theta == 0: # (12)
        return 0
    elif theta == 1: # (12)
        return 0
    else:
        return np.sin(theta * np.pi)


def iphi(z: float, omega1: int | float, omega2: int | float) -> float:
    theta = omega2 / abs(omega1)
    if z < 0:
        return - theta * phi(- z / theta)   # (10)
    elif z > 1:
        return - phi(1 + theta - theta * z) / theta     # (11)
    else:
        return phi(z)   # in [0,1] = \phi


def w_xt(u_lst: np.ndarray, x_lst: np.ndarray, t_lst: np.ndarray, omega: tuple):
    xi = 0
    d_omega = omega[1] - omega[0]
    for x in x_lst:
        ti = 0
        for t in t_lst:
            arg1 = x + t / omega[1]
            arg2 = x + t / omega[0]
            u_lst[ti][xi] = (omega[1] * phi(arg1) - omega[0] * phi(arg2)) / d_omega     # (9)
            ti += 1
        xi += 1
    return u_lst
