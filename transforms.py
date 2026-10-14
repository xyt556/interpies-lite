"""
interpies_lite.transforms
==========================

Pure-numpy/scipy implementations of potential-field transforms:
derivatives (Savitzky-Golay / finite-sample / Fourier), vertical derivative
via the ISVD (integrated second vertical derivative) method, upward
continuation, vertical integration, detrending and padding helpers.

No rasterio / GDAL / obspy dependency here.
"""
from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike
from scipy import ndimage as ndi
from scipy import signal
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures


def replace_edges(data: np.ndarray, ncells: int = 1) -> np.ndarray:
    out = data.copy()
    for i in range(ncells):
        out[i, :] = out[ncells, :]
        out[-1 - i, :] = out[-1 - ncells, :]
        out[:, i] = out[:, ncells]
        out[:, -1 - i] = out[:, -1 - ncells]
    return out


def fill_nodata(data: np.ndarray, invalid: np.ndarray | None = None) -> np.ndarray:
    if invalid is None:
        invalid = ~np.isfinite(data)
    if not invalid.any():
        return data.copy()
    return ndi.generic_filter(
        np.where(invalid, np.nan, data),
        lambda v: np.nanmean(v),
        size=3,
        mode="nearest",
    )


def simple_resample(data: np.ndarray, sampling: int = 2) -> np.ndarray:
    return data[::sampling, ::sampling]


def stats(data: np.ndarray) -> dict:
    arr = data[np.isfinite(data)]
    return {
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "median": float(np.median(arr)),
        "p02": float(np.percentile(arr, 2)),
        "p98": float(np.percentile(arr, 98)),
    }


def find_trend(point_coords: ArrayLike, data: ArrayLike, degree: int = 1, return_model: bool = False):
    X = np.asarray(point_coords)
    y = np.asarray(data).ravel()
    mask = np.isfinite(y)
    model = Pipeline([("poly", PolynomialFeatures(degree=degree)), ("reg", LinearRegression())])
    model.fit(X[mask], y[mask])
    if return_model:
        return model
    return model.predict(X).reshape(np.asarray(data).shape)


def savgol2d(degree: int, window_size: int) -> np.ndarray:
    if window_size % 2 == 0:
        raise ValueError("window_size must be odd")
    half = window_size // 2
    x = np.arange(-half, half + 1)
    order = np.arange(degree + 1)
    V = x[:, None] ** order[None, :]
    coef_0 = np.linalg.pinv(V)[0]
    kernel = np.outer(coef_0, coef_0)
    return kernel / kernel.sum()


def savgol_smooth(data: np.ndarray, deg: int = 3, win: int = 5, do_edges: bool = False) -> np.ndarray:
    kernel = savgol2d(deg, win)
    out = signal.convolve2d(np.nan_to_num(data), kernel, mode="same", boundary="symm")
    if do_edges:
        out = replace_edges(out, ncells=win // 2)
    return out


def savgol_deriv(
    data: np.ndarray,
    cellsize: float,
    direction: Literal["dx", "dy"] = "dx",
    deg: int = 3,
    win: int = 5,
    do_edges: bool = True,
) -> np.ndarray:
    if win % 2 == 0:
        raise ValueError("win must be odd")
    c_dx = signal.savgol_coeffs(win, deg, deriv=1, use="dot")
    c_smooth = signal.savgol_coeffs(win, deg, deriv=0, use="dot")
    if direction == "dx":
        kx, ky = c_dx.reshape(1, -1), c_smooth.reshape(-1, 1)
    else:
        kx, ky = c_smooth.reshape(1, -1), c_dx.reshape(-1, 1)
    kernel = ky @ kx
    out = signal.convolve2d(np.nan_to_num(data), kernel, mode="same", boundary="symm") / cellsize
    if do_edges:
        out = replace_edges(out, ncells=win // 2)
    return out


def fs_coefficients(tap: int = 5) -> np.ndarray:
    x = np.arange(-tap, tap + 1, dtype=float)
    A = np.vstack([np.ones_like(x), x]).T
    return np.linalg.pinv(A)[1]


def fs_deriv(data: np.ndarray, cellsize: float, direction: Literal["dx", "dy"] = "dx", tap: int = 5) -> np.ndarray:
    c = fs_coefficients(tap)
    kernel = c.reshape(1, -1) if direction == "dx" else c.reshape(-1, 1)
    return signal.convolve2d(np.nan_to_num(data), kernel, mode="same", boundary="symm") / cellsize


def next_pow2(x: int) -> int:
    p = 1
    while p < x:
        p <<= 1
    return p


def getk(nx: int, ny: int, dx: float, dy: float):
    fx = np.fft.fftfreq(nx, d=dx)
    fy = np.fft.fftfreq(ny, d=dy)
    kx, ky = np.meshgrid(2 * np.pi * fx, 2 * np.pi * fy)
    return kx, ky, np.sqrt(kx**2 + ky**2)


def pad_full(data: np.ndarray, mode: str = "reflect") -> np.ndarray:
    ny, nx = data.shape
    return np.pad(data, ((ny // 2, ny // 2), (nx // 2, nx // 2)), mode=mode)


def unpad_full(data: np.ndarray, nrows: int, ncols: int) -> np.ndarray:
    return data[nrows // 2:nrows // 2 + nrows, ncols // 2:ncols // 2 + ncols]


def fourier_transform(
    data: np.ndarray,
    cellsize: float,
    trans: Literal["dx", "dy", "dxdy", "dz", "vi", "upcont"] = "dx",
    order: int = 1,
    z: float = 0.0,
    eps: float = 1e-6,
    pad: bool = True,
) -> np.ndarray:
    work = pad_full(data) if pad else data.copy()
    ny, nx = work.shape
    kx, ky, k = getk(nx, ny, cellsize, cellsize)
    F = np.fft.fft2(np.nan_to_num(work))

    if trans == "dx":
        factor = (1j * kx) ** order
    elif trans == "dy":
        factor = (1j * ky) ** order
    elif trans == "dxdy":
        factor = (1j * kx) ** order * (1j * ky) ** order
    elif trans == "dz":
        factor = k**order
    elif trans == "vi":
        k_safe = k.copy()
        k_safe[k_safe == 0] = eps
        factor = k_safe ** (-order)
    elif trans == "upcont":
        factor = np.exp(-z * k)
    else:
        raise ValueError(f"Unknown transform: {trans}")

    out = np.real(np.fft.ifft2(F * factor))
    if pad:
        out = unpad_full(out, data.shape[0], data.shape[1])
    return out


def isvd(
    data: np.ndarray,
    cellsize: float,
    method: Literal["SG", "FS", "fourier"] = "SG",
    order: int = 1,
    deg: int = 4,
    win: int = 5,
    fs_tap: int = 5,
    do_edges: bool = True,
    **kwargs,
) -> np.ndarray:
    if method.upper() == "SG":
        dx2 = savgol_deriv(data, cellsize, "dx", deg=deg, win=win, do_edges=do_edges)
        dx2 = savgol_deriv(dx2, cellsize, "dx", deg=deg, win=win, do_edges=do_edges)
        dy2 = savgol_deriv(data, cellsize, "dy", deg=deg, win=win, do_edges=do_edges)
        dy2 = savgol_deriv(dy2, cellsize, "dy", deg=deg, win=win, do_edges=do_edges)
    elif method.upper() == "FS":
        dx2 = fs_deriv(fs_deriv(data, cellsize, "dx", tap=fs_tap), cellsize, "dx", tap=fs_tap)
        dy2 = fs_deriv(fs_deriv(data, cellsize, "dy", tap=fs_tap), cellsize, "dy", tap=fs_tap)
    else:
        dx2 = fourier_transform(data, cellsize, trans="dx", order=2, **kwargs)
        dy2 = fourier_transform(data, cellsize, trans="dy", order=2, **kwargs)

    dz2 = -(dx2 + dy2)
    if order == 2:
        return dz2
    if order == 1:
        return fourier_transform(dz2, cellsize, trans="vi", order=1, pad=True)
    raise NotImplementedError(f"ISVD order={order} not implemented")


def gauss(data: np.ndarray, sigma: float = 1.0) -> np.ndarray:
    return ndi.gaussian_filter(np.nan_to_num(data), sigma=sigma)


def laplacian(data: np.ndarray, cellsize: float) -> np.ndarray:
    return ndi.laplace(np.nan_to_num(data)) / (cellsize**2)
