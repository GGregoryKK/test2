def walsh_matrix(ca: list, n: int) -> list:
    """Compute the walsh matrix after `n` transformation with the given matrix.

    Args:
        ca: given square matrix.
        n: transformation count.

    Returns:
        Walsh matrix.
    """
    if ca:
        for _ in range(n):
            size = len(ca)
            res = [[0 for _ in range(size << 1)] for _ in range(size << 1)]
            for irow in range(size):
                for icol in range(size):
                    icol <<= 1
                    res[irow][icol] = ca[irow][icol >> 1]
                    res[irow][icol + 1] = ca[irow][icol >> 1]

                    res[size + irow][icol] = res[irow][icol] * (-1) ** icol
                    res[size + irow][icol + 1] = res[irow][icol + 1] * (-1) ** (icol + 1)
            ca = res
        return ca
    else:
        return ca


def fwf(sa: list) -> list:
    """Compute the one-dimensional discrete Walsh-Fourier Transform.

    Args:
        sa: one-dimensional list.

    Returns:
        Walsh-Fourier coefficients.
    """
    size = len(sa)
    step = size >> 1
    cp = sa.copy()

    if size == 2:
        return [(sa[0] + sa[1]) / 2, (sa[0] - sa[1]) / 2]
    else:
        for k in range(0, step):
            cp[k] = .5 * (sa[k * 2] + sa[k * 2 + 1])
            cp[step + k] = .5 * (sa[k * 2] - sa[k * 2 + 1])
        return fwf(cp[:step]) + fwf(cp[step:])


if __name__ == "__main__":
    arr: list = [[1]]
    rang_2: list = walsh_matrix(arr, 1)
    rang_4: list = walsh_matrix(arr, 2)
    print(f"\nrang 2:\n{rang_2}\n\nrang 4:\n{rang_4}\n\n")
