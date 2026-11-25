"""Inputs on a projected grid: bedrock topography (ETOPO 2022 with ice sheets removed, or free 10-30 m terrain
tiles for fine grids) and present climate (WorldClim worldwide, PRISM 800 m in the contiguous US)."""

import io
import time
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, map_coordinates

R = 6371000.0
LAPSE = 0.0065  # K per m, as in the ice model
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


ETOPO = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_bed_elev_netcdf/"
         "ETOPO_2022_v1_60s_N90W180_bed.nc.ascii")


def _get(url, tries=4):
    """Download with retries: the NOAA server sometimes drops long requests."""
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=900).read()
        except OSError:
            if attempt == tries - 1:
                raise
            time.sleep(10 * (attempt + 1))


def fetch_bed(lon_min, lon_max, lat_min, lat_max, stride):
    """Download (and cache) a box of ETOPO 2022 1' bedrock, every `stride` cells. Returns lon, lat, elevation."""
    i0, i1 = (int((v + 90) * 60) for v in (lat_min, lat_max))
    j0, j1 = (int((v + 180) * 60) for v in (lon_min, lon_max))
    i1, j1 = min(i1, 180 * 60 - 1), min(j1, 360 * 60 - 1)
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"bed_{i0}_{i1}_{j0}_{j1}_{stride}.npy"
    if not path.exists():
        print(f"downloading ETOPO 2022 bedrock -> {path.name}")
        text = _get(f"{ETOPO}?z%5B{i0}:{stride}:{i1}%5D%5B{j0}:{stride}:{j1}%5D").decode()
        rows = text.split("z.z[")[1].split("\n\n")[0].splitlines()[1:]
        np.save(path, np.array([[float(v) for v in r.split(",")[1:]] for r in rows]))
    lat = -90 + (np.arange(i0, i1 + 1, stride) + 0.5) / 60
    lon = -180 + (np.arange(j0, j1 + 1, stride) + 0.5) / 60
    return lon, lat, np.load(path)


TERRARIUM = "https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{}/{}/{}.png"


def fetch_terrain(lon_min, lon_max, lat_min, lat_max, z):
    """Mosaic of AWS Terrain Tiles (free open data, ~10-30 m sources) at zoom z, cached.
    Returns (first tile x, first tile y, elevation array). This is today's surface: use it only where there
    is no ice today, where it equals bedrock."""
    n = 2**z
    tx = lambda lon: int((lon + 180) / 360 * n)
    ty = lambda lat: int((1 - np.arcsinh(np.tan(np.radians(lat))) / np.pi) / 2 * n)
    x0, x1, y0, y1 = tx(lon_min), tx(lon_max), ty(lat_max), ty(lat_min)
    path = CACHE / f"terrain_{z}_{x0}_{x1}_{y0}_{y1}.npy"
    if not path.exists():
        print(f"downloading {(x1 - x0 + 1) * (y1 - y0 + 1)} terrain tiles (zoom {z})")
        def tile(xy):
            px = np.asarray(Image.open(io.BytesIO(urllib.request.urlopen(TERRARIUM.format(z, *xy), timeout=120).read())))
            return px[..., 0] * 256.0 + px[..., 1] + px[..., 2] / 256 - 32768
        with ThreadPoolExecutor(8) as pool:
            rows = [np.hstack(list(pool.map(tile, [(i, j) for i in range(x0, x1 + 1)]))) for j in range(y0, y1 + 1)]
        np.save(path, np.vstack(rows).astype(np.float32))
    return x0, y0, np.load(path)


def load_bed(proj, x, y, dx):
    """Bedrock elevation (m) on the projected grid defined by 1-D x, y arrays: ETOPO 2022 for ice-sheet grids,
    terrain tiles (pixels at most dx/2) for finer ones."""
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    if dx < 2000:
        z = int(np.ceil(np.log2(2 * 156543.03 * np.cos(np.radians(lat.mean())) / dx)))
        x0, y0, z_arr = fetch_terrain(lon.min() - 0.05, lon.max() + 0.05, lat.min() - 0.05, lat.max() + 0.05, z)
        col = ((lon + 180) / 360 * 2**z - x0) * 256 - 0.5
        row = ((1 - np.arcsinh(np.tan(np.radians(lat))) / np.pi) / 2 * 2**z - y0) * 256 - 0.5
        return map_coordinates(z_arr, [row, col], order=1, mode="nearest").astype(np.float64)
    elon, elat, z = fetch_bed(max(lon.min() - 0.5, -180), min(lon.max() + 0.5, 180),
                              max(lat.min() - 0.5, -90), min(lat.max() + 0.5, 90), max(1, int(dx / 1850 / 1.5)))
    rc = [(lat - elat[0]) / (elat[1] - elat[0]), (lon - elon[0]) / (elon[1] - elon[0])]
    return map_coordinates(z, rc, order=1, mode="nearest")


WORLDCLIM = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_{}_{}.zip"


def _worldclim(var, per_deg, rows, cols):
    """A window of the WorldClim 2.1 grids at 10 arc-minutes (per_deg 6) or 2.5 (24), oceans filled from the
    nearest land."""
    res = {6: "10m", 24: "2.5m"}[per_deg]
    path = CACHE / "worldclim" / (f"{var}.zip" if per_deg == 6 else f"{var}_{res}.zip")
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading WorldClim {var} ({res})")
        urllib.request.urlretrieve(WORLDCLIM.format(res, var), path)
    with zipfile.ZipFile(path) as z:
        a = np.stack([np.asarray(Image.open(io.BytesIO(z.read(n))))[rows, cols] for n in sorted(z.namelist()) if n.endswith(".tif")])
    _, (i, j) = distance_transform_edt(a[0] < -1e4, return_indices=True)  # no-data value differs per variable
    return a[:, i, j].astype(np.float64)


PRISM = "https://data.prism.oregonstate.edu/normals/us/800m/{0}/monthly/prism_{0}_us_30s_2020{1:02d}_avg_30y.zip"
PRISM_WEST, PRISM_NORTH = -125.0208333, 49.9375  # corner of the 30" grid (7025 x 3105)


def _prism(rows, cols):
    """PRISM 1991-2020 monthly normals (tmean °C, ppt mm), each (12, rows, cols) for a window of the grid; NaN off
    land. The whole grid is downloaded once (24 files of about 58 MB) and kept as two arrays on disk, 1 GB each."""
    out = []
    for var in ("tmean", "ppt"):
        path = CACHE / "prism" / f"{var}.npy"
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            print(f"downloading PRISM 800 m {var} normals (12 x ~58 MB)")
            def month(m):
                with zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(PRISM.format(var, m), timeout=600).read())) as z:
                    return np.asarray(Image.open(io.BytesIO(z.read(next(n for n in z.namelist() if n.endswith(".tif"))))))
            with ThreadPoolExecutor(6) as pool:
                months = list(pool.map(month, range(1, 13)))
            part = path.with_suffix(".part.npy")  # written under another name, so a half-written file is never used
            np.save(part, np.array(months, np.float32))
            part.rename(path)
        a = np.load(path, mmap_mode="r")[:, rows, cols].astype(np.float64)
        out.append(np.where(a < -9000, np.nan, a))
    return out


# PRISM (1991-2020) minus WorldClim (1970-2000) sea-level temperature by month, averaged over Idaho (K).
# Removed from PRISM so the LGM anomalies apply to the same baseline everywhere.
PRISM_SHIFT = np.array([1.04, 0.09, 0.53, 0.22, 0.36, 0.14, 1.11, 0.77, 0.72, -0.33, 0.03, 0.31])


def load_climate(proj, x, y, bed=None):
    """Present-day monthly temperature (°C), precipitation (mm) and their reference elevation.

    With `bed` (the bedrock on this grid, for mountain regions), WorldClim is read at 4.6 km and PRISM's 800 m
    normals replace it inside the contiguous US: they resolve mountain snowfall that coarser cells smear into
    the valleys.
    """
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    per_deg = 6 if bed is None else 24  # 18 km cells for the ice sheets, 4.6 km for mountain regions
    r, c = (90 - lat) * per_deg - 0.5, (lon + 180) * per_deg - 0.5
    rows, cols = (slice(max(int(v.min()) - 2, 0), int(v.max()) + 3) for v in (r, c))
    sample = lambda a, rc: np.array([map_coordinates(m, rc, order=1, mode="nearest", cval=np.nan) for m in a])
    T, P, z = (sample(_worldclim(v, per_deg, rows, cols), [r - rows.start, c - cols.start]) for v in ("tavg", "prec", "elev"))
    z_ref = z[0]
    if bed is not None:
        r, c = (PRISM_NORTH - lat) * 120 - 0.5, (lon - PRISM_WEST) * 120 - 0.5
        rows = slice(max(int(r.min()) - 2, 0), min(int(r.max()) + 3, 3105))
        cols = slice(max(int(c.min()) - 2, 0), min(int(c.max()) + 3, 7025))
        if rows.start < rows.stop and cols.start < cols.stop:
            pt, pp = (sample(a, [r - rows.start, c - cols.start]) for a in _prism(rows, cols))
            ok = np.isfinite(pt).all(0) & np.isfinite(pp).all(0)
            # PRISM is fitted to an 800 m DEM, which the bedrock matches on land.
            T = np.where(ok, pt - PRISM_SHIFT[:, None, None], T)
            P, z_ref = np.where(ok, pp, P), np.where(ok, bed, z_ref)
    return T, P, z_ref
