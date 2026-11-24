"""Model regions and their inputs: bedrock, LGM climate (with any fitted correction) and edge ice."""

from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from . import dem, export, pmip

WEB = Path(__file__).resolve().parent.parent / "web"

# name: (centre lon, centre lat, half-width km, half-height km, grid km, years)
# Mountain glaciers settle within 3,000 years; the ice-sheet lobes entering Idaho take far longer.
REGIONS = {
    "north_america": (-95.0, 58.0, 3300, 2600, 20, 40000),
    "idaho":         (-114.2, 45.6, 320, 440, 0.5, 15000),
    "lost_river":    (-113.65, 44.05, 45, 65, 0.2, 3000),
    "wind_river":    (-109.55, 42.95, 75, 95, 0.2, 3000),
    "sawtooth":      (-114.95, 44.10, 40, 50, 0.2, 3000),
}
# Nested regions take ice flowing in across their northern edge from a coarser run (Cordilleran lobes
# entering Idaho from British Columbia). Other edges stay ice-free: the coarse run over-glaciates the
# lowlands there.
PARENT = {"idaho": "north_america"}


def correction_path(region):
    return dem.CACHE / f"correction_{region}.npz"


def correction(region, lon, lat):
    """Fitted LGM temperature correction (K) at lon/lat (python -m glaciers.fit), or 0 if the region has none.
    The 2-3° PMIP4 models resolve the ice sheets but not mountain ranges, whose cooling they get wrong."""
    if not correction_path(region).exists():
        return 0.0
    c = np.load(correction_path(region))
    rc = [(lat - c["lats"][0]) / (c["lats"][1] - c["lats"][0]), (lon - c["lons"][0]) / (c["lons"][1] - c["lons"][0])]
    return map_coordinates(c["dT"], rc, order=1, mode="nearest")


def setup(region, use_correction=True):
    """Everything a run needs: grid, bedrock, LGM climate (T, P, their elevation) and edge ice."""
    lon0, lat0, hw, hh, km, years = REGIONS[region]
    dx = km * 1000.0
    proj = dem.LAEA(lon0, lat0)
    x = np.arange(-hw * 1000, hw * 1000 + 1, dx)
    y = np.arange(-hh * 1000, hh * 1000 + 1, dx)
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    bed0 = dem.load_bed(proj, x, y, dx)
    T, P, z_ref = dem.load_climate(proj, x, y, bed0 if dx < 2000 else None)  # PRISM for fine grids
    lgm_dT, ratio = pmip.sample(lon, lat)
    fitted = correction(region, lon, lat) if use_correction else 0.0
    climate = T + lgm_dT + fitted, P * ratio, z_ref  # present climate + LGM anomalies

    H_edge = None
    if region in PARENT:
        s, ice = export.read(WEB, PARENT[region], lon, lat)
        H_edge = np.where(ice, np.maximum(s - bed0, 0), 0)
        H_edge[:-1] = 0  # only the last row, the northern edge
    return dict(proj=proj, x=x, y=y, dx=dx, lon=lon, lat=lat, bed0=bed0, climate=climate, H_edge=H_edge,
                years=years, corrected=bool(np.any(fitted)))
