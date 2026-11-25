"""Model regions and their inputs: bedrock, LGM climate (with any fitted correction) and edge ice."""

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt, map_coordinates

from . import dem, export, outlines, pmip

WEB = Path(__file__).resolve().parent.parent / "web"

# name: (centre lon, centre lat, half-width km, half-height km, grid km, years)
# Mountain glaciers settle within 3,000 years; ice-sheet lobes take far longer.
REGIONS = {
    "north_america": (-95.0, 58.0, 3300, 2600, 20, 40000),
    "alaska":        (-151.0, 63.0, 850, 650, 1.5, 15000),
}
# The western United States is covered patch by patch: one 200 m region around each group of mapped mountain
# glaciers (python -m glaciers.west draws them up), rather than one grid over ground that was mostly ice-free.
PATCHES = Path(__file__).with_name("patches.json")
if PATCHES.exists():
    REGIONS.update({name: (*box, 0.2, 3000) for name, box in json.loads(PATCHES.read_text()).items()})
PATCH_NAMES = set(REGIONS) - {"north_america", "alaska"}
# In the patches the ice-sheet outline (Batchelor) counts as evidence only east of this longitude and south of this
# latitude, where it shows the mountain ice of central Idaho and Yellowstone. Further west it covers the Cascades
# whole; further north it is the Cordilleran ice sheet, which flowed in from Canada and which no patch can grow.
SHEETS = -120.0, 46.5
# Regions that take ice flowing in across their edges from a coarser run: Alaska is joined to the ice sheet.
PARENT = {"alaska": "north_america"}
# Mountains as they stood at the LGM: (summit lon, lat, height m, slope of the cone put back on top). Rainier lost
# its top to the Osceola collapse 5,600 years ago; before that it is thought to have reached about 4,900 m.
REBUILT = [(-121.7603, 46.8523, 4900.0, 0.7)]
# Where the outlines say nothing: (west, east, south, north, above this elevation m). In Washington's Cascades
# and Olympics they draw only the valley tongues (the share of ground mapped as ice falls with height there, where
# in Oregon and the Sierra it rises), and leave blank the high ground the tongues flowed from. Such ground is left
# out of fits and scores.
UNMAPPED = [(-125.0, -120.0, 45.8, 49.5, 1100.0)]

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
    """Everything a run needs, as a dict:

    proj, x, y, dx, lon, lat   the grid (projection, cell centres in metres, spacing, and their lon/lat)
    bed0, rebuilt              the ground, and how much of it was put back on top of today's (see REBUILT)
    climate, natural           ice-age (monthly temperature, monthly precipitation, their elevation); `natural`
                               is the same without mapped mode's changes, for colouring the ice (snow or bare)
    mapped, known              the field evidence on the grid, and where it can be compared with
    H_edge, years, corrected   ice thickness held at the edges, simulated years, whether a fit was applied

    In mapped mode the ice is tied to the field evidence instead of predicted: it is colder inside the mapped
    extent and nothing survives outside it, so ice fills the outlines and only its thickness and flow are
    computed.
    """
    # The grid: square cells on an equal-area projection centred on the region (x east, y north, metres).
    lon0, lat0, hw, hh, km, years = REGIONS[region]
    dx = km * 1000.0
    proj = dem.LAEA(lon0, lat0)
    x = np.arange(-hw * 1000, hw * 1000 + 1, dx)
    y = np.arange(-hh * 1000, hh * 1000 + 1, dx)
    X, Y = np.meshgrid(x, y)
    lon, lat = proj.inverse(X, Y)

    # The ground, and today's climate on it (fine grids get the finer climate data).
    bed0 = dem.load_bed(proj, x, y, dx)
    T, P, z_ref = dem.load_climate(proj, x, y, bed0 if dx < 2000 else None)
    away = lambda lon_p, lat_p: np.hypot(*np.subtract((X, Y), np.reshape(proj.forward(lon_p, lat_p), (2, 1, 1))))  # m
    known = np.ones(bed0.shape, bool)  # where the outlines can be compared with
    for west, east, south, north, above in UNMAPPED:
        known &= (lon < west) | (lon > east) | (lat < south) | (lat > north) | (bed0 < above)
    rebuilt = np.zeros_like(bed0)  # height put back on top of today's ground (m)
    for lon_p, lat_p, top, slope in REBUILT:  # a cone over the summit, wherever it stands above today's ground
        rebuilt = np.maximum(rebuilt, top - slope * away(lon_p, lat_p) - bed0)
    bed0 = bed0 + rebuilt

    # The ice-age climate: today's, plus the climate models' change, plus this region's fitted correction.
    lgm_dT, ratio = pmip.sample(lon, lat)
    fitted = correction(region, lon, lat) if use_correction else 0.0
    natural = climate = T + lgm_dT + fitted, P * ratio, z_ref

    # The field evidence on this grid and, in mapped mode, the climate bent to it: colder inside the mapped ice
    # (building up over MARGIN from its outline), far too warm for ice outside.
    # (For the patches, the mountain-glacier outlines only: the ice-sheet outline covers the Cascades whole.)
    evidence = outlines.rasterize(proj, x, y)
    if region in PATCH_NAMES:  # the western patches use the ice-sheet outline only where it shows mountain ice
        evidence = outlines.rasterize(proj, x, y, sheets=False) | evidence & (lon > SHEETS[0]) & (lat < SHEETS[1])
    if mapped:
        inside = np.minimum(distance_transform_edt(evidence) * dx / MARGIN, 1)
        climate = climate[0] + np.where(inside > 0, MAPPED[0] * inside, MAPPED[1]), *climate[1:]

    # Ice entering from a coarser run: its thickness along this region's edges (only the outermost cells are used).
    H_edge = None
    if region in PARENT:
        s, ice = export.read(WEB, PARENT[region], lon, lat)
        H_edge = np.where(ice, np.maximum(s - bed0, 0), 0)
    return dict(proj=proj, x=x, y=y, dx=dx, lon=lon, lat=lat, bed0=bed0, rebuilt=rebuilt, climate=climate,
                natural=natural, mapped=evidence, known=known, H_edge=H_edge, years=years,
                corrected=bool(np.any(fitted)))
