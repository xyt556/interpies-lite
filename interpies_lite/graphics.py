"""
interpies_lite.graphics
=======================

Plotting: hillshaded 2D grids with configurable colormap normalisation,
contours and colorbar.
"""
from __future__ import annotations

from typing import Literal

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, ListedColormap, BoundaryNorm

from interpies_lite import colors as icolors
from interpies_lite import transforms as T


# ---------------------------------------------------------------------------
# normalisation helpers
# ---------------------------------------------------------------------------

def _auto_norm(data: np.ndarray, mode: str = "equalize") -> Normalize:
    finite = data[np.isfinite(data)]
    if mode == "equalize":
        # histogram-equalise the lookup: resample colors at quantiles
        lo, hi = np.percentile(finite, [2, 98])
        return Normalize(vmin=lo, vmax=hi)
    if mode == "autolevels":
        lo, hi = np.percentile(finite, [2, 98])
        return Normalize(vmin=lo, vmax=hi)
    if mode == "linear":
        return Normalize(vmin=finite.min(), vmax=finite.max())
    raise ValueError(mode)


# ---------------------------------------------------------------------------
# hillshade
# ---------------------------------------------------------------------------

def _hillshade(
    data: np.ndarray,
    cellsize: float,
    zf: float = 10.0,
    azdeg: float = 45.0,
    altdeg: float = 45.0,
) -> np.ndarray:
    """Simple hillshade; returns values in [0, 1]."""
    az = np.radians(360.0 - azdeg)
    alt = np.radians(altdeg)
    dy, dx = np.gradient(np.nan_to_num(data) * zf, cellsize, cellsize)
    slope = np.pi / 2.0 - np.arctan(np.sqrt(dx**2 + dy**2))
    aspect = np.arctan2(-dx, dy)
    shaded = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect)
    return np.clip(shaded, 0, 1)


# ---------------------------------------------------------------------------
# main show function
# ---------------------------------------------------------------------------

def show_grid(
    grid,
    ax=None,
    cmap: str = "geosoft",
    cmap_norm: Literal["equalize", "autolevels", "linear"] = "equalize",
    hs: bool = True,
    zf: float = 10.0,
    azdeg: float = 45.0,
    altdeg: float = 45.0,
    hs_contrast: float = 1.5,
    cmap_brightness: float = 1.0,
    alpha: float = 0.7,
    contours: bool | int = False,
    colorbar: bool = True,
    cb_contours: bool = False,
    std_range: float = 1.0,
    figsize: tuple[float, float] = (8, 8),
    title: str | None = None,
    **kwargs,
):
    """Display a Grid with hillshading and optional contours.

    Parameters mirror the original interpies ``Grid.show()``.
    """
    data = grid.data
    extent = grid.extent
    cmap_obj = icolors.load_cmap(cmap)

    # optional brightness tweak
    if cmap_brightness != 1.0:
        arr = cmap_obj(np.linspace(0, 1, 256))
        arr[:, :3] = np.clip(arr[:, :3] * cmap_brightness, 0, 1)
        cmap_obj = ListedColormap(arr)

    norm = _auto_norm(data, cmap_norm)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # colour layer
    im = ax.imshow(
        np.ma.masked_invalid(data),
        extent=extent,
        origin="upper",
        cmap=cmap_obj,
        norm=norm,
        alpha=alpha if hs else 1.0,
        **kwargs,
    )

    # hillshade overlay
    if hs:
        shade = _hillshade(data, grid.cellsize, zf=zf, azdeg=azdeg, altdeg=altdeg)
        shade = np.clip((shade - 0.5) * hs_contrast + 0.5, 0, 1)
        ax.imshow(
            shade,
            extent=extent,
            origin="upper",
            cmap="gray",
            alpha=1.0 - alpha,
            zorder=im.get_zorder() - 0.1,
        )

    # contours
    if contours:
        n_levels = 32 if contours is True else int(contours)
        cs = ax.contour(
            np.ma.masked_invalid(data),
            extent=extent,
            origin="upper",
            levels=n_levels,
            colors="k",
            linewidths=0.4,
            alpha=0.6,
        )
        if cb_contours:
            # draw contour level ticks on the colorbar later
            pass

    if colorbar:
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(grid.name)

    ax.set_title(title or grid.name)
    ax.set_xlabel("Easting")
    ax.set_ylabel("Northing")
    return ax


def save_grid_image(grid, output_file: str, dpi: int = 150, **kwargs) -> None:
    """Render the grid without title/colorbar/annotations and save to file."""
    kwargs.setdefault("colorbar", False)
    kwargs.setdefault("title", "")
    ax = show_grid(grid, **kwargs)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.figure.savefig(output_file, dpi=dpi, bbox_inches="tight", pad_inches=0)
    plt.close(ax.figure)
