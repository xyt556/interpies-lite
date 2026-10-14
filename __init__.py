"""
interpies_lite
==============

A from-scratch, dependency-light re-implementation of the key features of
the original ``interpies`` library: loading geophysical grids, applying
potential-field transforms (derivatives, ISVD, upward continuation,
edge-enhancement filters like tilt / HGM / TDR / local wavenumber) and
plotting with hillshading.

Example
-------
>>> import interpies_lite as ip
>>> g = ip.open("data/brtpgrd.gxf")
>>> g.info()
>>> g.detrend().tilt().show()
"""
from __future__ import annotations

from interpies_lite.grid import Grid
from interpies_lite.io import open_rasterio as open

__version__ = "0.1.0"
__all__ = ["Grid", "open", "open_rasterio"]
