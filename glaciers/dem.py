"""Inputs on a projected grid: bedrock topography (ETOPO 2022, ice sheets removed) and present climate (WorldClim)."""

import io
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, map_coordinates

R = 6371000.0
CACHE = Path(__file__).resolve().parent.parent / "data"


class LAEA:
    """Spherical Lambert azimuthal equal-area projection centred on (lon0, lat0)."""

    def __init__(self, lon0, lat0):
        self.lon0, self.lat0 = np.radians(lon0), np.radians(lat0)

    def forward(self, lon, lat):
        lam, phi = np.radians(lon) - self.lon0, np.radians(lat)
        s0, c0 = np.sin(self.lat0), np.cos(self.lat0)
        k = np.sqrt(2 / (1 + s0 * np.sin(phi) + c0 * np.cos(phi) * np.cos(lam)))
        x = R * k * np.cos(phi) * np.sin(lam)
        y = R * k * (c0 * np.sin(phi) - s0 * np.cos(phi) * np.cos(lam))
        return x, y

    def inverse(self, x, y):
        rho = np.maximum(np.hypot(x, y), 1e-9)
        c = 2 * np.arcsin(np.clip(rho / (2 * R), -1, 1))
        s0, c0 = np.sin(self.lat0), np.cos(self.lat0)
        lat = np.arcsin(np.cos(c) * s0 + y * np.sin(c) * c0 / rho)
        lon = self.lon0 + np.arctan2(x * np.sin(c), rho * c0 * np.cos(c) - y * s0 * np.sin(c))
        return (np.degrees(lon) + 180) % 360 - 180, np.degrees(lat)


OPENDAP = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/"
           "60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc.ascii")


def fetch_bed(lon_min, lon_max, lat_min, lat_max, stride):
    """Download (and cache) an ETOPO 2022 1-arc-minute bedrock box. Returns lon, lat, elevation."""
    i0, i1 = (int((v + 90) * 60) for v in (lat_min, lat_max))
    j0, j1 = (int((v + 180) * 60) for v in (lon_min, lon_max))
    i1, j1 = min(i1, 10799), min(j1, 21599)
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"bed_{i0}_{i1}_{j0}_{j1}_{stride}.npy"
    if not path.exists():
        print(f"downloading ETOPO 2022 bedrock -> {path.name}")
        url = f"{OPENDAP}?z%5B{i0}:{stride}:{i1}%5D%5B{j0}:{stride}:{j1}%5D"
        text = urllib.request.urlopen(url, timeout=600).read().decode()
        rows = text.split("z.z[")[1].split("\n\n")[0].splitlines()[1:]
        np.save(path, np.array([[float(v) for v in r.split(",")[1:]] for r in rows]))
    z = np.load(path)
    lat = -90 + (np.arange(i0, i1 + 1, stride) + 0.5) / 60
    lon = -180 + (np.arange(j0, j1 + 1, stride) + 0.5) / 60
    return lon, lat, z


def load_bed(proj, x, y, dx):
    """Bedrock elevation (m) on the projected grid defined by 1-D x, y arrays."""
    X, Y = np.meshgrid(x, y)
    lon, lat = proj.inverse(X, Y)
    stride = max(1, int(dx / 1852 / 1.5))  # 1 arc-minute is ~1852 m
    pad = 0.5
    elon, elat, z = fetch_bed(
        max(lon.min() - pad, -180), min(lon.max() + pad, 180),
        max(lat.min() - pad, -90), min(lat.max() + pad, 90), stride,
    )
    rows = (lat - elat[0]) / (elat[1] - elat[0])
    cols = (lon - elon[0]) / (elon[1] - elon[0])
    return map_coordinates(z, [rows, cols], order=1, mode="nearest")


WORLDCLIM = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_{}.zip"


def _worldclim(var):
    """WorldClim 2.1 10-arc-minute grids (n, 1080, 2160) with oceans filled from the nearest land, and the ocean mask."""
    path = CACHE / "worldclim" / f"{var}.zip"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading WorldClim {var}")
        urllib.request.urlretrieve(WORLDCLIM.format(var), path)
    with zipfile.ZipFile(path) as z:
        a = np.stack([np.asarray(Image.open(io.BytesIO(z.read(n)))) for n in sorted(z.namelist()) if n.endswith(".tif")])
    ocean = a[0] < -1e4  # no-data value differs per variable
    _, (i, j) = distance_transform_edt(ocean, return_indices=True)
    return a[:, i, j].astype(np.float64), ocean


def load_climate(proj, x, y):
    """Present-day monthly temperature (°C), precipitation (mm), their reference elevation, and the ocean mask."""
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    rc = [(90 - lat) * 6 - 0.5, (lon + 180) * 6 - 0.5]
    sample = lambda a: map_coordinates(a, rc, order=1, mode="nearest")
    T, P = (np.array([sample(m) for m in _worldclim(v)[0]]) for v in ("tavg", "prec"))
    elev, ocean = _worldclim("elev")
    return T, P, sample(elev[0]), map_coordinates(ocean.astype(float), rc, order=0, mode="nearest") > 0.5
