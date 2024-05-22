import numpy as np
from numpy.fft import fft2, ifft2
from pypher.pypher import psf2otf


def l0(s, lambda_value=0.0015, betamax=1e5, kappa=1.05):
    s = np.array(s, dtype=float)
    fx = np.array([[1, -1]])
    fy = np.array([[1], [-1]])
    otfFx = psf2otf(fx, (N, M))
    otfFy = psf2otf(fy, (N, M))
    ns = fft2(S)
    sq_otf = np.abs(otfFx) ** 2 + np.abs(otfFy) ** 2
    beta = 2 * lambda_value
    while beta < betamax:
        denom = 1 + beta * sq_otf
        h = np.vstack((np.diff(s, 1, 0), (s[0, :] - s[-1, :])))
        v = np.hstack((np.diff(s, 1, 1), (s[:, 0] - s[:, -1]).reshape(-1, 1)))
        t = np.array(h ** 2 + v ** 2 < lambda_value / beta, dtype=bool)
        h[t], v[t] = 0, 0
        nomin = np.vstack(((h[-1, :] - h[0, :]), -np.diff(h, 1, 0))) + np.hstack(
            ((v[:, -1] - v[:, 0]).reshape(-1, 1), -np.diff(v, 1, 1)))
        FS = (ns + beta * fft2(nomin)) / denom
        s = ifft2(FS).real
        beta = beta * kappa
    return S
