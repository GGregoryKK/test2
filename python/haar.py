import numpy as np
from matplotlib import pyplot as plt

from pywt import wavedec


def conv(f, g):
    result_size = len(f)
    if result_size <= 1:
        result_size = 2
    f_fft = np.fft.fft(f.copy(), result_size)
    g_fft = np.fft.fft(g.copy(), result_size)
    result_fft = f_fft * g_fft
    result = np.fft.ifft(result_fft).real
    return result


def rec(data, h=None, g=None, lvl=1):
    hsize = len(data) // 2 ** lvl
    data[:hsize * 2] = conv(np.insert(data[:hsize], np.arange(1, hsize + 1), 0), h) + \
                       conv(np.insert(data[hsize: hsize * 2], np.arange(1, hsize + 1), 0), g)
    for l in range(1, lvl):
        hsize *= 2
        data[:hsize * 2] = conv(np.insert(data[:hsize], np.arange(1, hsize + 1), 0), h) + \
                           conv(np.insert(data[hsize: hsize * 2], np.arange(1, hsize + 1), 0), g)
    return data


def dec(v, lvl=1, h=None, g=None):
    size = len(v) >> 1
    a = conv(v, h)[1::2]
    v = np.concatenate((a, conv(v, g)[1::2]), axis=None)
    for l in range(1, lvl):
        a = conv(a, h)[1::2]
        v = np.concatenate((a, conv(v[:size], g)[1::2], v[size:]), axis=None)
        size >>= 1
    return v


def fht(arr):
    step = len(arr) >> 1
    cp = arr.copy()
    c = 1 / np.sqrt(2)
    for _ in range(int(np.log2(len(arr)))):
        for k in range(0, step):
            cp[k] = c * (arr[k * 2] + arr[k * 2 + 1])
            cp[step + k] = c * (arr[k * 2] - arr[k * 2 + 1])
        step >>= 1
        arr = cp.copy()
    return cp


def ifht(arr):
    step = 1
    cp = arr.copy()
    c = 1 / np.sqrt(2)
    for _ in range(int(np.log2(len(arr)))):
        for k in range(0, step):
            cp[k * 2] = c * (arr[k] + arr[step + k])
            cp[k * 2 + 1] = c * (arr[k] - arr[step + k])
        step <<= 1
        arr = cp.copy()
    return arr


def to_string(a):
    return ', '.join([str(round(i, 4)) for i in a])


def present_haar_transform():
    N = 8
    v = np.arange(N, dtype=float)

    z = 1 / np.sqrt(2)
    h = [z, z]
    g = [-z, z]
    print(f"Изначальные данные:\n\t{to_string(v)}")
    print(f"Прямое БПХ:\n\t{to_string(fht(v))}")
    print(f"Обратное БПХ:\n\t{to_string(ifht(fht(v)))}")

    print(f"Прямое ПХ чз свертку:\n\t{to_string(dec(v, lvl=3, h=h, g=g))}")
    print(f"Обратное ПХ чз свертку:\n\t{to_string(rec(fht(v), lvl=3, h=h, g=g[::-1]))}")
    print(f"Готовая реализация:\n\t{wavedec(v, 'haar')}")


def graph_ifht(x, arr, plot_steps=None):
    step = 1
    cp = arr.copy()
    c = 1 / np.sqrt(2)
    for _ in range(int(np.log2(len(arr)))):
        a = []
        for k in range(0, step):
            cp[k * 2] = c * (arr[k] + arr[step + k])
            cp[k * 2 + 1] = c * (arr[k] - arr[step + k])
            a.append(c * (arr[k] + arr[step + k]))
            a.append(c * (arr[k] - arr[step + k]))
        step <<= 1
        if step in plot_steps:
            plt.step(x, np.repeat(np.array(a, dtype=float), len(x) // len(a)), where='mid', linestyle='-.',
                     label="L " + str(step))
        arr = cp.copy()
    return arr


def wavelet_series():
    N = 2 ** 3
    x = np.linspace(-10, 10, N)
    func = lambda y: y ** 3
    v = func(x)
    plt.plot(x, v, label="Signal")
    plt.plot(x, ifht(fht(v.copy())), label="Recon", linestyle='--')
    graph_ifht(x, fht(v.copy()), plot_steps=[2, 4, 8, 16])
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    present_haar_transform()
    wavelet_series()
