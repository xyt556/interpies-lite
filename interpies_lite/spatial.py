"""
interpies_lite.spatial
======================

Lightweight coordinate helpers.  Uses pyproj when available; falls back to
a no-op with a warning.
"""
from __future__ import annotations

import numpy as np


def project_points(points, src_srs: str = "EPSG:4326", dst_srs: str = "EPSG:32633"):
    """Reproject an (N, 2) array of points from src_srs to dst_srs."""
    try:
        from pyproj import Transformer
    except ImportError:
        raise RuntimeError("pyproj is required: pip install pyproj")
    t = Transformer.from_crs(src_srs, dst_srs, always_xy=True)
    pts = np.asarray(points, dtype=float)
    x, y = t.transform(pts[:, 0], pts[:, 1])
    return np.column_stack([x, y])


def grid_extent(west: float, north: float, cellsize: float, nrows: int, ncols: float):
    """Return matplotlib-style [left, right, bottom, top]."""
    return [west, west + ncols * cellsize, north - nrows * cellsize, north]
