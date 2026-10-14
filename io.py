"""
interpies_lite.io
=================

Thin wrapper around rasterio for opening single-band rasters.
"""
from __future__ import annotations

from interpies_lite.grid import Grid


def open_rasterio(path: str, band: int = 1, name: str | None = None) -> Grid:
    try:
        import rasterio
    except ImportError as e:
        raise RuntimeError(
            "rasterio is required to open rasters. Install with `pip install rasterio`."
        ) from e
    with rasterio.open(path) as ds:
        return Grid.from_dataset(ds, band=band, name=name)
