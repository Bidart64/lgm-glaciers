"""Model regions and their inputs: bedrock, LGM climate (with any fitted correction) and edge ice."""

from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates

from . import dem, export, pmip

WEB = Path(__file__).resolve().parent.parent / "web"

# name: (centre lon, centre lat, half-width km, half-height km, grid km, years, time step yr, extra cooling K)
# The 2-3° PMIP4 models resolve the ice sheets but not mountain ranges, whose LGM cooling they underestimate.
# Idaho and the Wind River Range instead get corrections fitted to the mapped ice (python -m glaciers.fit <region>).
REGIONS = {
    "alps":          (9.5, 46.6, 600, 380, 2, 3000, 10, -3.5),
    "pyrenees":      (0.5, 42.6, 300, 150, 2, 3000, 10, -3.5),
    "new_zealand":   (171.0, -43.0, 400, 450, 3, 3000, 10, -3.5),
    "patagonia":     (-71.0, -47.0, 450, 1100, 4, 5000, 10, -3.5),
    "europe":        (20.0, 65.0, 2800, 2200, 20, 40000, 25, 0),
    "north_america": (-95.0, 58.0, 3300, 2600, 20, 40000, 25, 0),
    "idaho":         (-114.2, 45.6, 320, 440, 0.5, 3000, 10, 0),
    "lost_river":    (-113.65, 44.05, 45, 65, 0.2, 3000, 5, 0),
    "wind_river":    (-109.55, 42.95, 75, 95, 0.2, 3000, 5, 0),
}
# Regions that reuse another region's fitted correction (it is stored in lon/lat).
CORRECTION_FROM = {"lost_river": "idaho"}
# Nested regions take ice flowing in across one edge from a coarser run (Cordilleran lobes entering Idaho
# from British Columbia). Other edges stay ice-free: the coarse run over-glaciates the lowlands there.
PARENT = {"idaho": ("north_america", "north")}


def correction_path(region):
    return dem.CACHE / f"correction_{region}.npz"


def correction(region, lon, lat):
    """Fitted LGM temperature correction (K) at lon/lat, or 0 if the region has none."""
    path = correction_path(CORRECTION_FROM.get(region, region))
    if not path.exists():
        return 0.0
    c = np.load(path)
    lons, lats = c["lons"], c["lats"]
    rc = [(lat - lats[0]) / (lats[1] - lats[0]), (lon - lons[0]) / (lons[1] - lons[0])]
    return map_coordinates(c["dT"], rc, order=1, mode="nearest")


def setup(region, km=None, dT=None, use_correction=True):
    """Everything a run needs: grid, bedrock, LGM climate and edge ice. `km` overrides the grid spacing."""
    lon0, lat0, hw, hh, region_km, years, dt, region_dT = REGIONS[region]
    dT = region_dT if dT is None else dT
    dx = (km or region_km) * 1000.0
    proj = dem.LAEA(lon0, lat0)
    x = np.arange(-hw * 1000, hw * 1000 + 1, dx)
    y = np.arange(-hh * 1000, hh * 1000 + 1, dx)
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    bed0 = dem.load_bed(proj, x, y, dx)
    T, P, z_ref, ocean = dem.load_climate(proj, x, y, bed0 if dx < 2000 else None)  # PRISM for fine grids
    lgm_dT, ratio = pmip.sample(lon, lat)
    fitted = correction(region, lon, lat) if use_correction else 0.0
    climate = T + lgm_dT + dT + fitted, P * ratio, z_ref, ocean  # present climate + LGM anomalies

    H_edge = None
    if region in PARENT:
        parent, side = PARENT[region]
        s, ice = export.read(WEB, parent, lon, lat)
        H_edge = np.where(ice, np.maximum(s - bed0, 0), 0)
        if side == "north":
            H_edge[:-1] = 0  # last row is the northern edge
    return dict(proj=proj, x=x, y=y, dx=dx, lon=lon, lat=lat, bed0=bed0, climate=climate, H_edge=H_edge,
                years=years, dt=dt, dT=dT, corrected=bool(np.any(fitted)))
