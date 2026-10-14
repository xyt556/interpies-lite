"""
interpies_lite.colors
=====================

Minimal built-in colormaps: Geosoft clra-style, parula, and access to
matplotlib built-ins.  Avoids shipping the 35 KB lookup tables of the
original library; perceptually uniform maps from matplotlib (viridis,
plasma, inferno, magma, cividis) are recommended.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Geosoft "clra" approximate ramp: blue -> cyan -> green -> yellow -> red
_CLRA_NODES = [
    (0.00, (0.00, 0.00, 0.50)),
    (0.15, (0.00, 0.00, 1.00)),
    (0.35, (0.00, 0.75, 1.00)),
    (0.50, (0.00, 1.00, 0.50)),
    (0.65, (0.75, 1.00, 0.00)),
    (0.85, (1.00, 0.75, 0.00)),
    (1.00, (0.50, 0.00, 0.00)),
]

# MATLAB parula key stops
_PARULA_NODES = [
    (0.00, (0.208, 0.166, 0.529)),
    (0.25, (0.00, 0.387, 0.882)),
    (0.50, (0.00, 0.68, 0.6)),
    (0.75, (0.65, 0.82, 0.25)),
    (1.00, (0.98, 0.90, 0.13)),
]


def _make(name: str, nodes) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list(name, nodes)


# register builtins
_BUILTINS = {
    "geosoft": _make("geosoft", _CLRA_NODES),
    "clra": _make("clra", _CLRA_NODES),
    "parula": _make("parula", _PARULA_NODES),
}


def load_cmap(name: str = "geosoft"):
    """Return a colormap by name; falls back to matplotlib built-ins."""
    if name in _BUILTINS:
        return _BUILTINS[name]
    if name in plt.colormaps():
        return plt.get_cmap(name)
    raise ValueError(f"Unknown colormap: {name}")


def list_cmaps() -> list[str]:
    return sorted(list(_BUILTINS.keys()) + [c for c in plt.colormaps() if not c.endswith("_r")])
