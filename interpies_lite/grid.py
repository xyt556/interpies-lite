"""
interpies_lite.grid
===================

The :class:`Grid` object: a 2D array with georeferencing, plus all the
potential-field methods (derivatives, enhancements, continuation) that make
interpies useful.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from interpies_lite import transforms as T


@dataclass
class Grid:
    """A single-band georeferenced raster.

    Attributes
    ----------
    data : 2D float array (nodata -> NaN)
    transform : 6-tuple or affine.Affine
        Affine mapping (row, col) -> (x, y).
    crs : str or None
        Coordinate reference system (EPSG:xxxx string or PROJ.4).
    name, filename : str
    nodata : float or None
    """

    data: np.ndarray
    transform: Any = None
    crs: str | None = None
    name: str = "Unknown"
    filename: str = "Unknown"
    nodata: float | None = None

    # ---------- post-init derived attributes ----------
    def __post_init__(self) -> None:
        self.data = np.asarray(self.data, dtype=float)
        if self.nodata is not None:
            self.data[self.data == self.nodata] = np.nan
        if self.transform is None:
            # default: origin (0,0), 100 m cells
            self.transform = (100.0, 0.0, 0.0, 0.0, -100.0, 0.0)

    # ---------- derived properties ----------
    @property
    def nrows(self) -> int:
        return self.data.shape[0]

    @property
    def ncols(self) -> int:
        return self.data.shape[1]

    @property
    def cellsize(self) -> float:
        return abs(float(self.transform[0]))

    @property
    def y_cellsize(self) -> float:
        return abs(float(self.transform[4]))

    @property
    def x_origin(self) -> float:
        """X of upper-left corner."""
        return float(self.transform[2])

    @property
    def y_origin(self) -> float:
        """Y of upper-left corner."""
        return float(self.transform[5])

    @property
    def extent(self) -> list[float]:
        """matplotlib-style extent: [left, right, bottom, top]."""
        left = self.x_origin
        right = left + self.ncols * self.cellsize
        top = self.y_origin
        bottom = top - self.nrows * self.y_cellsize
        return [left, right, bottom, top]

    @property
    def mask(self) -> np.ndarray:
        return ~np.isfinite(self.data)

    # ---------- constructors ----------
    @classmethod
    def from_array(
        cls,
        array: np.ndarray,
        west: float = 0.0,
        north: float = 0.0,
        cellsize: float = 100.0,
        y_cellsize: float | None = None,
        crs: str | None = None,
        name: str = "Unknown",
    ) -> "Grid":
        ys = y_cellsize if y_cellsize is not None else cellsize
        transform = (cellsize, 0.0, west, 0.0, -ys, north)
        return cls(data=np.asarray(array), transform=transform, crs=crs, name=name)

    @classmethod
    def from_dataset(cls, dataset, band: int = 1, name: str | None = None) -> "Grid":
        """Build a Grid from an open rasterio dataset."""
        data = dataset.read(band).astype(float)
        crs = dataset.crs.to_string() if dataset.crs is not None else None
        return cls(
            data=data,
            transform=dataset.transform,
            crs=crs,
            name=name or os.path.splitext(os.path.basename(dataset.name))[0],
            filename=dataset.name,
            nodata=dataset.nodata,
        )

    # ---------- IO ----------
    def save(self, output_file: str) -> None:
        """Write as GeoTIFF."""
        try:
            import rasterio
        except ImportError as e:
            raise RuntimeError(
                "rasterio is required for save(); install with `pip install rasterio`"
            ) from e
        profile = {
            "driver": "GTiff",
            "height": self.nrows,
            "width": self.ncols,
            "count": 1,
            "dtype": "float64",
            "crs": self.crs,
            "transform": self.transform,
        }
        with rasterio.open(output_file, "w", **profile) as dst:
            dst.write(np.nan_to_num(self.data), 1)

    # ---------- geometry / preprocessing ----------
    def _new(self, data: np.ndarray, name: str | None = None, **kw) -> "Grid":
        base = dict(transform=self.transform, crs=self.crs, filename=self.filename, nodata=self.nodata)
        base.update(kw)
        return Grid(data=data, name=name or self.name, **base)

    def clip(self, xmin: float, xmax: float, ymin: float, ymax: float) -> "Grid":
        left, right, bottom, top = self.extent
        c1 = max(0, int(np.floor((xmin - left) / self.cellsize)))
        c2 = min(self.ncols, int(np.ceil((xmax - left) / self.cellsize)))
        r1 = max(0, int(np.floor((top - ymax) / self.y_cellsize)))
        r2 = min(self.nrows, int(np.ceil((top - ymin) / self.y_cellsize)))
        sub = self.data[r1:r2, c1:c2]
        new_left = left + c1 * self.cellsize
        new_top = top - r1 * self.y_cellsize
        new_transform = (self.transform[0], 0.0, new_left, 0.0, self.transform[4], new_top)
        return self._new(sub, name=f"{self.name}_clip", transform=new_transform)

    def resample(self, sampling: int = 2) -> "Grid":
        sub = T.simple_resample(self.data, sampling)
        new_cs = self.cellsize * sampling
        new_ycs = self.y_cellsize * sampling
        new_transform = (new_cs, 0.0, self.x_origin, 0.0, -new_ycs, self.y_origin)
        return self._new(sub, name=f"{self.name}_res{sampling}", transform=new_transform)

    def scale(self, factor: float) -> "Grid":
        return self._new(self.data * factor, name=f"{self.name}_scaled")

    def detrend(self, degree: int = 1, sampling: int = 4) -> "Grid":
        sub = T.simple_resample(np.nan_to_num(self.data), sampling)
        ny, nx = sub.shape
        ys = np.linspace(self.y_origin, self.y_origin - ny * self.y_cellsize * sampling, ny)
        xs = np.linspace(self.x_origin, self.x_origin + nx * self.cellsize * sampling, nx)
        XX, YY = np.meshgrid(xs, ys)
        trend_sub = T.find_trend(np.c_[XX.ravel(), YY.ravel()], sub.ravel(), degree=degree)
        # resize trend back to original shape
        from scipy.ndimage import zoom
        zy = self.nrows / ny
        zx = self.ncols / nx
        trend_full = zoom(trend_sub.reshape(ny, nx), (zy, zx), order=1)
        return self._new(self.data - trend_full, name=f"{self.name}_detrend{degree}")

    def fill_nodata(self) -> "Grid":
        return self._new(T.fill_nodata(self.data), name=f"{self.name}_filled")

    def smooth(self, method: str = "SG", deg: int = 3, win: int = 5, sigma: float = 1.0) -> "Grid":
        if method.upper() == "SG":
            out = T.savgol_smooth(self.data, deg=deg, win=win)
        else:
            out = T.gauss(self.data, sigma=sigma)
        return self._new(out, name=f"{self.name}_smooth")

    def laplacian(self) -> "Grid":
        return self._new(T.laplacian(self.data, self.cellsize), name=f"{self.name}_lap")

    def info(self) -> None:
        s = T.stats(self.data)
        print(f"Grid: {self.name}")
        print(f"  file      : {self.filename}")
        print(f"  shape     : {self.nrows} rows x {self.ncols} cols")
        print(f"  cellsize  : {self.cellsize} x {self.y_cellsize}")
        print(f"  extent    : {self.extent}")
        print(f"  crs       : {self.crs}")
        print(f"  min/max   : {s['min']:.3f} / {s['max']:.3f}")
        print(f"  mean/std  : {s['mean']:.3f} +/- {s['std']:.3f}")

    # ---------- derivatives ----------
    def dx(self, method: str = "SG", deg: int = 3, win: int = 5, fs_tap: int = 5) -> "Grid":
        if method.upper() == "FS":
            out = T.fs_deriv(self.data, self.cellsize, "dx", tap=fs_tap)
        elif method.lower() == "fourier":
            out = T.fourier_transform(self.data, self.cellsize, trans="dx", order=1)
        else:
            out = T.savgol_deriv(self.data, self.cellsize, "dx", deg=deg, win=win)
        return self._new(out, name=f"{self.name}_dx")

    def dy(self, method: str = "SG", deg: int = 3, win: int = 5, fs_tap: int = 5) -> "Grid":
        if method.upper() == "FS":
            out = T.fs_deriv(self.data, self.cellsize, "dy", tap=fs_tap)
        elif method.lower() == "fourier":
            out = T.fourier_transform(self.data, self.cellsize, trans="dy", order=1)
        else:
            out = T.savgol_deriv(self.data, self.cellsize, "dy", deg=deg, win=win)
        return self._new(out, name=f"{self.name}_dy")

    def dx2(self, method: str = "SG", deg: int = 4, win: int = 5, fs_tap: int = 5) -> "Grid":
        return self.dx(method=method, deg=deg, win=win, fs_tap=fs_tap).dx(method=method, deg=deg, win=win, fs_tap=fs_tap)

    def dy2(self, method: str = "SG", deg: int = 4, win: int = 5, fs_tap: int = 5) -> "Grid":
        return self.dy(method=method, deg=deg, win=win, fs_tap=fs_tap).dy(method=method, deg=deg, win=win, fs_tap=fs_tap)

    def dxdy(self, method: str = "SG", deg: int = 3, win: int = 5, fs_tap: int = 5) -> "Grid":
        return self.dx(method=method, deg=deg, win=win, fs_tap=fs_tap).dy(method=method, deg=deg, win=win, fs_tap=fs_tap)

    def dz(self, method: str = "fourier", order: int = 1, **kw) -> "Grid":
        if method.lower() == "fourier":
            out = T.fourier_transform(self.data, self.cellsize, trans="dz", order=order, **kw)
        else:
            out = T.isvd(self.data, self.cellsize, method=method, order=order, **kw)
        return self._new(out, name=f"{self.name}_dz{order}")

    def isvd(self, method: str = "SG", order: int = 1, deg: int = 3, win: int = 5, fs_tap: int = 5) -> "Grid":
        out = T.isvd(self.data, self.cellsize, method=method, order=order, deg=deg, win=win, fs_tap=fs_tap)
        return self._new(out, name=f"{self.name}_dz{order}")

    def vi(self, order: int = 1, eps: float = 1e-6) -> "Grid":
        out = T.fourier_transform(self.data, self.cellsize, trans="vi", order=order, eps=eps)
        return self._new(out, name=f"{self.name}_vi{order}")

    # ---------- enhancement transforms ----------
    def hgm(self, method: str = "SG", deg: int = 3, win: int = 5, fs_tap: int = 5) -> "Grid":
        gx = self.dx(method=method, deg=deg, win=win, fs_tap=fs_tap).data
        gy = self.dy(method=method, deg=deg, win=win, fs_tap=fs_tap).data
        return self._new(np.sqrt(gx * gx + gy * gy), name=f"{self.name}_HGM")

    def tilt(self, hgm_method: str = "SG", dz_method: str = "isvd",
             deg: int = 3, win: int = 5, fs_tap: int = 5, alpha: float = 1.0) -> "Grid":
        hgm_grid = self.hgm(method=hgm_method, deg=deg, win=win, fs_tap=fs_tap)
        if dz_method.lower() == "isvd":
            dz_grid = self.isvd(method=hgm_method, order=1, deg=deg, win=win, fs_tap=fs_tap)
        else:
            dz_grid = self.dz(method="fourier", order=1)
        out = np.arctan(alpha * dz_grid.data / hgm_grid.data) * 180.0 / np.pi
        return self._new(out, name=f"{self.name}_tilt")

    def tg(self, method: str = "SG", dz_method: str = "isvd",
           deg: int = 3, win: int = 5, fs_tap: int = 5) -> "Grid":
        gx = self.dx(method=method, deg=deg, win=win, fs_tap=fs_tap).data
        gy = self.dy(method=method, deg=deg, win=win, fs_tap=fs_tap).data
        if dz_method.lower() == "isvd":
            gz = self.isvd(method=method, order=1, deg=deg, win=win, fs_tap=fs_tap).data
        else:
            gz = self.dz(method="fourier", order=1).data
        return self._new(np.sqrt(gx**2 + gy**2 + gz**2), name=f"{self.name}_TG")

    def lw(self, method: str = "SG", dz_method: str = "isvd",
           deg: int = 4, win: int = 5, fs_tap: int = 5) -> "Grid":
        gx = self.dx(method=method, deg=deg, win=win, fs_tap=fs_tap).data
        gy = self.dy(method=method, deg=deg, win=win, fs_tap=fs_tap).data
        if dz_method.lower() == "isvd":
            gz = self.isvd(method=method, order=1, deg=deg, win=win, fs_tap=fs_tap).data
            gzz = self.isvd(method=method, order=2, deg=deg, win=win, fs_tap=fs_tap).data
        else:
            gz = self.dz(method="fourier", order=1).data
            gzz = self.dz(method="fourier", order=2).data
        # d/dx(dz) via a follow-up dx on gz
        gzx = T.savgol_deriv(gz, self.cellsize, "dx", deg=deg, win=win)
        gzy = T.savgol_deriv(gz, self.cellsize, "dy", deg=deg, win=win)
        num = gzx * gx + gzy * gy + gzz * gz
        den = gx * gx + gy * gy + gz * gz
        return self._new(num / den, name=f"{self.name}_LW")

    def tahg(self, method: str = "SG", dz_method: str = "isvd",
             deg: int = 4, win: int = 5, fs_tap: int = 5, alpha: float = 1.0) -> "Grid":
        hgm_grid = self.hgm(method=method, deg=deg, win=win, fs_tap=fs_tap)
        out = hgm_grid.tilt(hgm_method=method, dz_method=dz_method,
                           deg=deg, win=win, fs_tap=fs_tap, alpha=alpha)
        out.name = f"{self.name}_TAHG"
        return out

    # ---------- continuation / filtering ----------
    def up(self, z: float = 500.0) -> "Grid":
        out = T.fourier_transform(self.data, self.cellsize, trans="upcont", z=z)
        return self._new(out, name=f"{self.name}_UC{z}")

    def hp_filter_uc(self, z: float = 5000.0) -> "Grid":
        regional = T.fourier_transform(self.data, self.cellsize, trans="upcont", z=z)
        return self._new(self.data - regional, name=f"{self.name}_HPUC{z}")

    # ---------- display ----------
    def show(self, **kwargs):
        from interpies_lite.graphics import show_grid
        return show_grid(self, **kwargs)

    def save_image(self, output_file: str, **kwargs) -> None:
        from interpies_lite.graphics import save_grid_image
        save_grid_image(self, output_file, **kwargs)
