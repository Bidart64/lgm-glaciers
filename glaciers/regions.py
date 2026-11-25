"""Model regions and their inputs: bedrock, LGM climate (with any fitted correction) and edge ice."""

from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

from . import dem, export, outlines, pmip

WEB = Path(__file__).resolve().parent.parent / "web"

# name: (centre lon, centre lat, half-width km, half-height km, grid km, years)
# Mountain glaciers settle within 3,000 years; the ice-sheet lobes entering Idaho take far longer.
REGIONS = {
    "north_america": (-95.0, 58.0, 3300, 2600, 20, 40000),
    "idaho":         (-114.2, 45.6, 320, 440, 0.5, 15000),
    "lost_river":    (-113.65, 44.05, 45, 65, 0.2, 3000),
    "wind_river":    (-109.55, 42.95, 75, 95, 0.2, 3000),
    "sawtooth":      (-114.95, 44.10, 40, 50, 0.2, 3000),
    "rainier":       (-121.80, 46.80, 60, 50, 0.2, 3000),
    "pyrenees":      (0.60, 42.70, 220, 80, 0.3, 3000),
}
# Nested regions take ice flowing in across their northern edge from a coarser run (Cordilleran lobes
# entering Idaho from British Columbia). Other edges stay ice-free: the coarse run over-glaciates the
# lowlands there.
PARENT = {"idaho": "north_america"}
# Mountains as they stood at the LGM: (summit lon, lat, height m, slope of the cone put back on top). Rainier lost
# its top to the Osceola collapse 5,600 years ago; before that it is thought to have reached about 4,900 m.
REBUILT = {"rainier": (-121.7603, 46.8523, 4900.0, 0.7)}
# Regions scored against the mountain-glacier outlines only: the ice-sheet outline covers the Cascades whole.
NO_SHEETS = {"rainier"}
# Elevation (m) above which the outlines say nothing: around Rainier they draw only the valley tongues, and leave
# blank both the mountain itself, which still carries glaciers, and the slopes above the ice-age snowline.
# Ground above it is left out of fits and scores.
UNMAPPED_ABOVE = {"rainier": 1100.0}

MAPPED = -4.0, 20.0   # K added inside and outside the mapped ice in mapped mode: ice fills the outlines (more
                      # cooling adds little) and none survives outside
MARGIN = 1500.0       # inside, the cooling builds up over this distance from the outline, so snouts thin out (m)


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


def setup(region, use_correction=True, mapped=False):
    """Everything a run needs: grid, bedrock, LGM climate (T, P, their elevation) and edge ice.
    In mapped mode the ice is tied to the field evidence instead of predicted: it is colder inside the
    mapped extent and nothing survives outside it, so ice fills the outlines and only its thickness and flow
    are computed. `natural` is the climate without that, for colouring the ice (snow or bare)."""
    lon0, lat0, hw, hh, km, years = REGIONS[region]
    dx = km * 1000.0
    proj = dem.LAEA(lon0, lat0)
    x = np.arange(-hw * 1000, hw * 1000 + 1, dx)
    y = np.arange(-hh * 1000, hh * 1000 + 1, dx)
    lon, lat = proj.inverse(*np.meshgrid(x, y))
    bed0 = dem.load_bed(proj, x, y, dx)
    T, P, z_ref = dem.load_climate(proj, x, y, bed0 if dx < 2000 else None)  # PRISM for fine grids
    known = bed0 < UNMAPPED_ABOVE.get(region, np.inf)  # where the outlines can be compared with
    rebuilt = np.zeros_like(bed0)  # height put back on top of today's ground (m)
    if region in REBUILT:
        lon_s, lat_s, top, slope = REBUILT[region]
        rebuilt = np.maximum(top - slope * np.hypot(*np.subtract(np.meshgrid(x, y), np.reshape(proj.forward(lon_s, lat_s), (2, 1, 1)))) - bed0, 0)
        bed0 = bed0 + rebuilt
    lgm_dT, ratio = pmip.sample(lon, lat)
    fitted = correction(region, lon, lat) if use_correction else 0.0
    natural = climate = T + lgm_dT + fitted, P * ratio, z_ref  # present climate + LGM anomalies
    evidence = outlines.rasterize(proj, x, y, region not in NO_SHEETS)  # mapped LGM ice on this grid
    if mapped:
        inside = np.minimum(distance_transform_edt(evidence) * dx / MARGIN, 1)
        climate = climate[0] + np.where(inside > 0, MAPPED[0] * inside, MAPPED[1]), *climate[1:]

    H_edge = None
    if region in PARENT:
        s, ice = export.read(WEB, PARENT[region], lon, lat)
        H_edge = np.where(ice, np.maximum(s - bed0, 0), 0)
        H_edge[:-1] = 0  # only the last row, the northern edge
    return dict(proj=proj, x=x, y=y, dx=dx, lon=lon, lat=lat, bed0=bed0, rebuilt=rebuilt, climate=climate, natural=natural, mapped=evidence, known=known,
                H_edge=H_edge,
                years=years, corrected=bool(np.any(fitted)))
