"""Inputs on a projected grid: bedrock topography (ETOPO 2022, ice sheets removed) and present climate
(WorldClim worldwide, PRISM 800 m in the contiguous US)."""

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


OPENDAP = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/{s}s/"
           "{s}s_bed_elev_netcdf/ETOPO_2022_v1_{s}s_N90W180_bed.nc.ascii")
# 15" comes in 15° tiles named by their north-west corner. Its bedrock variant only covers the Arctic, so the
# surface variant is used: identical to bedrock wherever there is no ice sheet today.
OPENDAP_15S = ("https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/15s/15s_surface_elev_netcdf/"
               "ETOPO_2022_v1_15s_{}{:02d}{}{:03d}_surface.nc.ascii")


def _get(url, tries=4):
    """Download with retries: the NOAA server sometimes drops long requests."""
    for attempt in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=900).read()
        except OSError:
            if attempt == tries - 1:
                raise
            time.sleep(10 * (attempt + 1))


def _tiles(per_deg):
    """(first global row, first global column, rows, columns, url) of each file of an ETOPO 2022 resolution."""
    if per_deg < 240:
        yield 0, 0, 180 * per_deg, 360 * per_deg, OPENDAP.format(s=3600 // per_deg)
        return
    for north in range(-75, 91, 15):
        for west in range(-180, 180, 15):
            url = OPENDAP_15S.format("N" if north >= 0 else "S", abs(north), "W" if west < 0 else "E", abs(west))
            yield (north - 15 + 90) * 240, (west + 180) * 240, 3600, 3600, url


def fetch_bed(lon_min, lon_max, lat_min, lat_max, stride, per_deg=60):
    """Download (and cache) an ETOPO 2022 bedrock box at 1' (per_deg=60), 30" (120) or 15" (240).
    Returns lon, lat, elevation."""
    i0, i1 = (int((v + 90) * per_deg) for v in (lat_min, lat_max))
    j0, j1 = (int((v + 180) * per_deg) for v in (lon_min, lon_max))
    i1, j1 = min(i1, 180 * per_deg - 1), min(j1, 360 * per_deg - 1)
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"bed_{i0}_{i1}_{j0}_{j1}_{stride}{'' if per_deg == 60 else f'_{3600 // per_deg}s'}.npy"
    if not path.exists():
        print(f"downloading ETOPO 2022 bedrock -> {path.name}")
        z = np.full((len(range(i0, i1 + 1, stride)), len(range(j0, j1 + 1, stride))), np.nan)
        for ti, tj, n, m, url in _tiles(per_deg):
            # global indices on the stride lattice that fall inside this tile
            a, b = max(i0, ti + (i0 - ti) % stride), min(i1, ti + n - 1)
            c, d = max(j0, tj + (j0 - tj) % stride), min(j1, tj + m - 1)
            if a > b or c > d:
                continue
            q = f"?z%5B{a - ti}:{stride}:{b - ti}%5D%5B{c - tj}:{stride}:{d - tj}%5D"
            text = _get(url + q).decode()
            rows = text.split("z.z[")[1].split("\n\n")[0].splitlines()[1:]
            block = np.array([[float(v) for v in r.split(",")[1:]] for r in rows])
            z[(a - i0) // stride:(a - i0) // stride + len(block), (c - j0) // stride:(c - j0) // stride + block.shape[1]] = block
        np.save(path, z)
    z = np.load(path)
    lat = -90 + (np.arange(i0, i1 + 1, stride) + 0.5) / per_deg
    lon = -180 + (np.arange(j0, j1 + 1, stride) + 0.5) / per_deg
    return lon, lat, z


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
    """Bedrock elevation (m) on the projected grid defined by 1-D x, y arrays."""
    X, Y = np.meshgrid(x, y)
    lon, lat = proj.inverse(X, Y)
    if dx < 450:  # finer than ETOPO's best (15"): free terrain tiles with pixels at most dx/2
        z = int(np.ceil(np.log2(2 * 156543.03 * np.cos(np.radians(lat.mean())) / dx)))
        x0, y0, z_arr = fetch_terrain(lon.min() - 0.05, lon.max() + 0.05, lat.min() - 0.05, lat.max() + 0.05, z)
        col = ((lon + 180) / 360 * 2**z - x0) * 256 - 0.5
        row = ((1 - np.arcsinh(np.tan(np.radians(lat))) / np.pi) / 2 * 2**z - y0) * 256 - 0.5
        return map_coordinates(z_arr, [row, col], order=1, mode="nearest").astype(np.float64)
    per_deg = 240 if dx < 1000 else 120 if dx < 2000 else 60  # finer bedrock for finer grids
    stride = max(1, int(dx / (111e3 / per_deg) / 1.5))
    pad = 0.5
    elon, elat, z = fetch_bed(
        max(lon.min() - pad, -180), min(lon.max() + pad, 180),
        max(lat.min() - pad, -90), min(lat.max() + pad, 90), stride, per_deg,
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


PRISM = "https://data.prism.oregonstate.edu/normals/us/800m/{0}/monthly/prism_{0}_us_30s_2020{1:02d}_avg_30y.zip"
PRISM_WEST, PRISM_NORTH = -125.0208333, 49.9375  # corner of the 30" grid (7025 x 3105)


def _prism(rows, cols):
    """PRISM 1991-2020 monthly normals (tmean °C, ppt mm), each (12, rows, cols) for a window of the grid; NaN off land.
    Each month is about 58 MB; only the window is cached."""
    path = CACHE / f"prism_{rows.start}_{rows.stop}_{cols.start}_{cols.stop}.npz"
    if not path.exists():
        def month(args):
            var, m = args
            with zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(PRISM.format(var, m), timeout=600).read())) as z:
                a = np.asarray(Image.open(io.BytesIO(z.read(next(n for n in z.namelist() if n.endswith(".tif"))))))
            return np.where(a[rows, cols] < -9000, np.nan, a[rows, cols])
        print("downloading PRISM 800 m climate normals (24 x ~58 MB)")
        with ThreadPoolExecutor(6) as pool:
            a = np.array(list(pool.map(month, [(v, m) for v in ("tmean", "ppt") for m in range(1, 13)])))
        np.savez_compressed(path, T=a[:12], P=a[12:])
    d = np.load(path)
    return d["T"], d["P"]


def load_climate(proj, x, y, bed=None, shift=None):
    """Present-day monthly temperature (°C), precipitation (mm), their reference elevation, and the ocean mask.

    With `bed` (the 30" bedrock on this grid), PRISM's 800 m normals replace WorldClim inside the contiguous US.
    They resolve mountain snowfall that WorldClim's 18 km cells smear into the valleys.
    Also returns the monthly PRISM-minus-WorldClim shift that was removed (None without PRISM); pass `shift`
    to reuse one, e.g. the shift a borrowed climate correction was fitted with.
    """
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    rc = [(90 - lat) * 6 - 0.5, (lon + 180) * 6 - 0.5]
    sample = lambda a: map_coordinates(a, rc, order=1, mode="nearest")
    T, P = (np.array([sample(m) for m in _worldclim(v)[0]]) for v in ("tavg", "prec"))
    elev, ocean = _worldclim("elev")
    z_ref, ocean = sample(elev[0]), map_coordinates(ocean.astype(float), rc, order=0, mode="nearest") > 0.5
    if bed is not None:
        r, c = (PRISM_NORTH - lat) * 120 - 0.5, (lon - PRISM_WEST) * 120 - 0.5
        rows = slice(max(int(r.min()) - 2, 0), min(int(r.max()) + 3, 3105))
        cols = slice(max(int(c.min()) - 2, 0), min(int(c.max()) + 3, 7025))
        if rows.start < rows.stop and cols.start < cols.stop:
            PT, PP = _prism(rows, cols)
            rc = [r - rows.start, c - cols.start]
            pt = np.array([map_coordinates(m, rc, order=1, mode="nearest", cval=np.nan) for m in PT])
            pp = np.array([map_coordinates(m, rc, order=1, mode="nearest", cval=np.nan) for m in PP])
            ok = np.isfinite(pt).all(0) & np.isfinite(pp).all(0)
            # PRISM is fitted to an 800 m DEM; the 30" bedrock matches it on land.
            # Shift its 1991-2020 temperatures to WorldClim's 1970-2000 baseline, month by month (regional
            # sea-level means), so the LGM anomalies and calibrations apply unchanged.
            if shift is None:
                shift = ((pt + LAPSE * bed) - (T + LAPSE * z_ref))[:, ok].mean(1)
            print(f"PRISM minus WorldClim baseline: {np.mean(shift):+.2f} K (removed)")
            T = np.where(ok, pt - shift[:, None, None], T)
            P, z_ref = np.where(ok, pp, P), np.where(ok, bed, z_ref)
    return T, P, z_ref, ocean, shift
