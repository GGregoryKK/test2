import numpy as np
import cmath as cm
from cmath import log, pi, exp


def bit_reverse(k, r):
    """Bit-reversed version of an integer number.

    Args:
        k: The number to be bit-reversed.
        r: The number of bits to take into consideration when reversing.

    Returns:
        The number k, bit-reversed according to integers with r bits.
    """
    l = 0
    for _ in np.arange(0, r):
        l = (l << 1) + (k & 1)
        k = (k >> 1)
    return l


def fft(sx, inv: bool = False):
    """Compute the one-dimensional discrete Fourier Transform.
        Args:
            sx: one-dimensional list.
            inv: True if FFT, Else if iFFT
        Returns:
            Permuted Fourier coefficients
        """
    size = len(sx)
    p2 = int(log(size, 2).real)
    c0 = 1 if not inv else 2
    c1 = -2j if not inv else 2j
    delta = size >> 1
    x = sx.copy()
    for _ in np.arange(p2):
        for b in np.arange(delta):
            for a in np.arange(0, size, delta * 2):
                j = a + b
                z = exp(c1 * pi * j / (delta << 1))
                t0 = (x[j] + x[j + delta]) / c0
                t1 = (x[j] - x[j + delta]) * z / c0
                x[j] = t0
                x[j + delta] = t1
        delta = delta >> 1
    for k in np.arange(0, size):
        l = bit_reverse(k, p2)
        sx[l] = x[k]
    return sx


def round_im(num: complex, ndig: int = 1) -> complex:
    """Round a complex number to a given precision in decimal digits.

    Args:
        num: given complex number.
        ndig: precision in decimal digits.

    Returns:
        rounded complex number.
    """
    return round(num.real, ndig) + round(num.imag, ndig) * 1j


def fft_rec(sa: list) -> list:
    """Compute the one-dimensional discrete Walsh-Fourier Transform.

    Args:
        sa: one-dimensional list.

    Returns:
        Permuted Fourier coefficients
    """
    size = len(sa)
    step = size >> 1
    cp = sa.copy()

    if size == 2:
        z = cm.exp((-2j * cm.pi) / step)
        return [round_im((sa[0] + sa[1])), round_im((sa[0] - sa[1]) * z)]

    else:
        for k in range(step):
            z = cm.exp((-1j * cm.pi * k) / step)
            cp[k] = round_im(sa[k] + sa[step + k])
            cp[step + k] = round_im(sa[k] * z - sa[step + k] * z)

        return fft_rec(cp[:step]) + fft_rec(cp[step:])
