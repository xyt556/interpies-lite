"""Smoke test: build a synthetic grid and run every transform end-to-end."""
import sys
from pathlib import Path

import numpy as np

# make the package importable when running from the repo root
sys.path.insert(0, str(Path(__file__).resolve().parent))

from interpies_lite import Grid
import interpies_lite as ip


def make_synthetic(n=120):
    x = np.linspace(-5000, 5000, n)
    y = np.linspace(-5000, 5000, n)
    XX, YY = np.meshgrid(x, y)
    # two Gaussian anomalies + a linear trend + noise
    g1 = 80 * np.exp(-((XX - 1500) ** 2 + (YY - 800) ** 2) / (2 * 800**2))
    g2 = -50 * np.exp(-((XX + 1800) ** 2 + (YY + 1200) ** 2) / (2 * 1000**2))
    trend = 0.002 * XX + 0.001 * YY
    rng = np.random.default_rng(0)
    data = g1 + g2 + trend + rng.normal(0, 0.5, size=XX.shape)
    return Grid.from_array(data, west=0, north=10000, cellsize=100, crs="EPSG:32633", name="synth")


def main():
    g = make_synthetic()
    print("== grid.info() ==")
    g.info()

    print("\n== derivatives ==")
    for m in [g.dx(), g.dy(), g.dx2(), g.dy2(), g.dxdy(), g.dz(method="fourier"), g.isvd()]:
        print(f"  {m.name:20s}  range=({np.nanmin(m.data):8.2f}, {np.nanmax(m.data):8.2f})")

    print("\n== enhancement transforms ==")
    for m in [g.hgm(), g.tilt(), g.tg(), g.lw(), g.tahg()]:
        print(f"  {m.name:20s}  range=({np.nanmin(m.data):8.2f}, {np.nanmax(m.data):8.2f})")

    print("\n== continuation / filtering ==")
    for m in [g.up(z=500), g.hp_filter_uc(z=2000), g.smooth(), g.laplacian(), g.vi()]:
        print(f"  {m.name:20s}  range=({np.nanmin(m.data):8.2f}, {np.nanmax(m.data):8.2f})")

    print("\n== preprocessing ==")
    g2 = g.detrend(degree=1).clip(1000, 8000, -8000, -1000).resample(2)
    print(f"  clipped/resampled shape: {g2.data.shape}, extent={g2.extent}")

    print("\n== plotting (headless) ==")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, m in zip(axes.ravel(), [g, g.hgm(), g.tilt(), g.tg(), g.up(500), g.hp_filter_uc(2000)]):
        m.show(ax=ax, hs=True, colorbar=False)
    out = Path(__file__).parent / "smoke_test_output.png"
    fig.savefig(out, dpi=100, bbox_inches="tight")
    print(f"  saved -> {out}")
    plt.close(fig)

    print("\nALL OK")


if __name__ == "__main__":
    main()
