import scipy
import numpy as np
import skimage

from image import GrayImage as GI
from assadasd import otf2psf, get_psf
from for_QA import add_gaussian_noise
import cv2
from pypher.pypher import zero_pad

eps = 1e-16


def psf2otf(psf, shape):
    if np.all(psf == 0):
        return np.zeros_like(psf)

    pad = False

    if len(shape) == 3:
        pad = True
        orig_shape = shape
        shape = (shape[0], shape[1])

    inshape = psf.shape
    psf = zero_pad(psf, shape, position='corner')

    for axis, axis_size in enumerate(inshape):
        psf = np.roll(psf, -int(axis_size / 2), axis=axis)

    if pad:
        psf2 = np.zeros(orig_shape)
        psf2[:, :, 0] = psf
    else:
        psf2 = psf
    otf = scipy.fft.fftn(psf2)

    n_ops = np.sum(psf.size * np.log2(psf.shape))
    otf = np.real_if_close(otf, tol=n_ops)

    return otf


def reflect(I):
    if len(I.shape) == 3:
        Is = I[::-1, ::-1, :]
    else:
        Is = I[::-1, ::-1]
    return Is


def estimate_noise(I):
    H = I.shape[0]
    W = I.shape[1]
    M = np.array([[1, -2, 1], [-2, 4, -2], [1, -2, 1]])
    S = np.sum(np.sum(np.abs(scipy.signal.convolve2d(I, M))))
    S = S * np.sqrt(0.5 * np.pi) / (6.0 * (W - 2.0) * (H - 2.0))
    return S


def estimate_nsr(I: np.ndarray):
    I = skimage.color.rgb2gray(I) if len(I.shape) == 3 else I
    en = estimate_noise(I)
    nsr = en ** 2 / np.var(I[:])
    return nsr


def aL(F, y, maxiter=500):
    L = y.copy()
    for i in range(maxiter):
        hp = y - F(L)
        hp = reflect(hp)
        d = (F(L + hp) - F(L - hp)) / 2.0
        d = reflect(d)
        L = L + d
    return L


def mLM(F, y, maxiter=500):
    nsr = estimate_nsr(y)
    a = 100.0 * nsr
    lm = y.copy()
    for i in range(maxiter):
        Flm = F(lm)
        H = scipy.fft.fftn(Flm) / (scipy.fft.fftn(lm) + eps)
        Hconj = H.conj
        num = Hconj * (scipy.fft.fftn(y - Flm))
        denom = (Hconj * H + a)
        lm = lm + np.real(scipy.fft.ifftn(num / denom))
    return lm


def mRL(F, y, maxiter=500):
    RL = y.copy()
    for i in range(maxiter):
        r1 = F(RL) / (np.abs(y) + eps)
        r1 = reflect(r1)
        r2 = F(r1)
        r2 = reflect(r2)
        RL = RL / (np.abs(r2) + eps)
        r1 = y / (np.abs(F(RL)) + eps)
        r1 = reflect(r1)
        r2 = F(r1)
        r2 = reflect(r2)
        RL = RL * r2
    return RL


def mW(F, y, maxiter=500):
    nsr = estimate_nsr(y)
    a = nsr
    W = y.copy()
    FW = F(W)
    H = scipy.fft.fftn(FW) / (scipy.fft.fftn(W) + eps)
    for i in range(1, maxiter + 1):
        H = H * (i - 1) / i + scipy.fft.fftn(FW) / (scipy.fft.fftn(W) + eps) / i
        Hconj = H.conj
        W = np.real(scipy.fft.ifftn(Hconj / (Hconj * H + a) * scipy.fft.fftn(y)))
        FW = F(W)
    return W


def pcVC(F, y, maxiter=500):
    TM = y.copy()
    for i in range(maxiter):
        H = scipy.fft.fftn(F(TM)) / (scipy.fft.fftn(TM) + eps)
        TM = TM + np.real(scipy.fft.ifftn((scipy.fft.fftn(y) / (H + eps) - scipy.fft.fftn(TM)) * np.absolute(H)))
    return TM


def pcVC_nsr(F, y, maxiter=500):
    nsr = estimate_nsr(y)
    a = 100.0 * nsr
    TM = y.copy()
    for i in range(maxiter):
        H = scipy.fft.fftn(F(TM)) / (scipy.fft.fftn(TM) + eps)
        Hconj = H.conj
        TM = TM + np.real(scipy.fft.ifftn(Hconj / (np.absolute(H) + a) * (scipy.fft.fftn(y) - H * scipy.fft.fftn(TM))))
    return TM
