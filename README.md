# interpies_lite

对原 `interpies` 库的轻量重实现。目标：解决原库依赖钉死（`gdal==3.8.4`）、久未更新的问题，保留核心功能与 API 习惯。

## 与原库的差异

| 维度 | 原 interpies | interpies_lite |
| --- | --- | --- |
| 栅格读写 | rasterio + 钉死 gdal==3.8.4 | rasterio>=1.3（可选，未装时只能用 `Grid.from_array`） |
| SEGY 模块 | obspy 依赖 | 砍掉 |
| 色表数据 | 35 KB 硬编码 LUT | 内置 Geosoft clra / parula 近似色表，其他直接用 matplotlib 色表 |
| Grid 实现 | 普通类 + 手写属性 | `@dataclass` + property，类型注解完整 |
| 依赖 | numpy/scipy/matplotlib/rasterio/gdal/scikit-image/scikit-learn | numpy/scipy/matplotlib/scikit-learn（rasterio/pyproj 可选） |
| Python | >=3.10 | >=3.10 |

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

## 安装

本包暂未发布到 PyPI，可通过以下任一方式安装：

### 方式 1：从 GitHub 直接安装（推荐）

公开仓库，无需下载，pip 自动拉取最新代码：

```bash
pip install git+https://github.com/xyt556/interpies-lite.git
```

指定版本/分支：

```bash
pip install git+https://github.com/xyt556/interpies-lite.git@main
```

### 方式 2：从本地目录安装

```bash
cd path/to/interpies_lite
pip install .
```

### 方式 3：从 wheel 文件安装

```bash
pip install dist/interpies_lite-0.1.0-py3-none-any.whl
```

### 方式 4：开发模式（可编辑，改代码即时生效）

```bash
cd path/to/interpies_lite
pip install -e .
```

### 可选依赖（读写 GeoTIFF/GXF）

```bash
pip install "interpies-lite[geo]"
# 等价于
pip install rasterio pyproj
```

核心依赖（必装）：numpy、scipy、matplotlib、scikit-learn。

## 冒烟测试

```bash
python smoke_test.py
# 在合成网格上跑通全部变换，输出 smoke_test_output.png
```
