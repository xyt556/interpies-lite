# interpies_lite

对原 `interpies` 库的轻量重实现。目标：解决原库依赖钉死（`gdal==3.8.4`）、久未更新的问题，保留核心功能与 API 习惯。

[![PyPI version](https://img.shields.io/pypi/v/interpies-lite.svg)](https://pypi.org/project/interpies-lite/)
[![Python versions](https://img.shields.io/pypi/pyversions/interpies-lite.svg)](https://pypi.org/project/interpies-lite/)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/xyt556/interpies-lite/blob/main/LICENSE)

## 与原库的差异

| 维度 | 原 interpies | interpies_lite |
| --- | --- | --- |
| 栅格读写 | rasterio + 钉死 gdal==3.8.4 | rasterio>=1.3（可选，未装时只能用 `Grid.from_array`） |
| SEGY 模块 | obspy 依赖 | 砍掉 |
| 色表数据 | 35 KB 硬编码 LUT | 内置 Geosoft clra / parula 近似色表，其他直接用 matplotlib 色表 |
| Grid 实现 | 普通类 + 手写属性 | `@dataclass` + property，类型注解完整 |
| 依赖 | numpy/scipy/matplotlib/rasterio/gdal/scikit-image/sklearn | numpy/scipy/matplotlib/sklearn（rasterio/pyproj 可选） |
| Python | >=3.10 | >=3.10 |

## 安装

**从 PyPI（推荐）：**

```bash
pip install interpies-lite
```

需要读写 GXF/GeoTIFF 时：

```bash
pip install "interpies-lite[geo]"
# 等价于
pip install rasterio pyproj
```

**从 GitHub 安装最新开发版：**

```bash
pip install git+https://github.com/xyt556/interpies-lite.git
```

**本地开发（可编辑模式）：**

```bash
git clone https://github.com/xyt556/interpies-lite.git
cd interpies-lite
pip install -e .[dev]
```

核心依赖（必装）：numpy、scipy、matplotlib、scikit-learn。

## 快速开始

```python
import interpies_lite as ip

# 从文件加载（需要 rasterio）
g = ip.open("data/brtpgrd.gxf")
g.info()

# 从 numpy 数组构造（无任何重依赖）
import numpy as np
data = np.random.randn(100, 100)
g = ip.Grid.from_array(data, west=0, north=0, cellsize=100, crs="EPSG:32633")

# 链式调用：去趋势 → 平滑 → 边界识别
out = g.detrend(degree=1).smooth(win=5).tilt()
out.show()

# 可用变换（全部返回新 Grid）：
#   dx / dy / dx2 / dy2 / dxdy / dz / isvd / vi
#   hgm / tilt / tg / lw / tahg
#   up(z) / hp_filter_uc(z) / laplacian
#   clip / resample / scale / detrend / fill_nodata / smooth
```

## 模块

```
interpies_lite/
├── __init__.py     # open() / Grid
├── grid.py         # Grid 类：IO、几何、预处理、全部位场方法
├── transforms.py   # 纯算法：导数、ISVD、傅里叶变换、延拓、平滑
├── graphics.py     # show / save_image：晕渲、色图、等值线
├── colors.py       # 内置色表
├── spatial.py      # pyproj 点投影（可选）
└── io.py           # rasterio 封装
```

## 示例

仓库自带两个 notebook：

| Notebook | 数据 | 内容 |
| --- | --- | --- |
| `notebooks/interpies_lite_demo.ipynb` | 合成数据 | 快速上手，无需 rasterio，开箱即跑 |
| `notebooks/interpies_lite_brtpgrd_demo.ipynb` | `data/brtpgrd.gxf`（真实航磁） | 课堂版完整教程，含算法原理与公式 |

`data/brtpgrd.gxf` 为 Colorado 地区航磁数据（707×629，50 m，NAD27/UTM13N），来自原 interpies 项目示例。

## 冒烟测试

```bash
python smoke_test.py
# 在合成网格上跑通全部变换，输出 smoke_test_output.png
```

## 发布

本项目通过 GitHub Actions 自动发布到 PyPI。打 tag 即触发：

```bash
git tag -a v0.1.0 -m "Release v0.1.0"
git push origin v0.1.0
```

## 许可证

BSD-3-Clause
