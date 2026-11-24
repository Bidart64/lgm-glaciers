"""Export a region's ice as one small PNG + metadata; the browser drapes it onto 10-30 m terrain."""

import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates

from .dem import LAEA
from .model import GLACIER

CREVASSE_STRAIN = 0.02  # 1/yr: tensile strain rate at which crevasses are fully developed


def flow_lines(vx, vy):
    """A value (0-1) that is constant along each flow line and varies smoothly across them, so its contours
    are flow lines: every cell takes smooth noise from where its ice came from, 256 cells upstream or the
    ice divide. Where tributaries join, their bands run side by side, like medial moraines."""
    speed = np.hypot(vx, vy)
    step = np.where(speed > 1e-3, 1 / np.maximum(speed, 1e-3), 0)
    pos = np.indices(vx.shape) - np.stack([vy, vx]) * step  # one cell upstream (rows run north)
    for _ in range(8):  # follow the path twice as far each pass
        pos = np.stack([map_coordinates(p, pos, order=1, mode="nearest") for p in pos])
    noise = gaussian_filter(np.random.default_rng(0).random(vx.shape), 4)
    noise = (noise - noise.min()) / np.ptp(noise)
    return map_coordinates(noise, pos, order=1, mode="nearest")


def write(web: Path, name, proj, x, y, dx, bed0, H, vx, vy, M, stats):
    """One PNG of three bands stacked north to south, each the size of the grid:
    1. R,G = ice surface + 10000 (m, 16-bit), B = glacier thickness (m, capped at 255; 0 under GLACIER)
    2. R = crevassing (largest tensile strain rate), G = direction of that stretching (0-255 for 0-180°,
       east = 0, anticlockwise; crevasses open across it), B = 128 + 25 × mass balance (m/yr): above 128 is
       snow-covered accumulation area, below is bare ice
    3. R = flow-line value (see flow_lines)

    The surface is relative to today's ground (no isostatic depression), so the browser can compare it
    directly with present-day elevation tiles. Thick ice is drawn up to that surface, flat across a valley;
    thin ice is draped over the fine terrain at the model's thickness.
    """
    ice = H > GLACIER
    dist, near = distance_transform_edt(~ice, return_indices=True)
    # Just outside the ice (valley walls), carry the nearest ice surface rather than the rock, so the
    # browser's interpolated surface stays flat until it meets the walls instead of ramping up them.
    surface = np.where(dist > 2, bed0, (bed0 + H)[tuple(near)])
    surface = np.clip(np.round(surface + 10000), 0, 65535).astype(np.uint32)

    dvx_dy, dvx_dx = np.gradient(vx, dx)
    dvy_dy, dvy_dx = np.gradient(vy, dx)
    shear, half = 0.5 * (dvx_dy + dvy_dx), 0.5 * (dvx_dx - dvy_dy)
    stretch = 0.5 * (dvx_dx + dvy_dy) + np.hypot(half, shear)  # largest principal strain rate
    angle = 0.5 * np.arctan2(shear, half) % np.pi               # its direction
    zero = np.zeros_like(H)
    bands = [[surface >> 8, surface & 255, np.where(ice, np.minimum(H, 255), 0)],
             [ice * np.clip(stretch / CREVASSE_STRAIN, 0, 1) * 255, angle / np.pi * 255 % 256, np.clip(128 + 25 * M, 0, 255)],
             [flow_lines(vx, vy) * 255, zero, zero]]
    px = np.vstack([np.round(np.stack(b, -1))[::-1] for b in bands]).astype(np.uint8)  # row 0 = north
    out = web / "data"
    out.mkdir(parents=True, exist_ok=True)
    Image.fromarray(px).save(out / f"{name}.png", optimize=True)

    lon, lat = proj.inverse(*np.meshgrid(x, y))
    meta = {"name": name, "lon0": float(np.degrees(proj.lon0)), "lat0": float(np.degrees(proj.lat0)),
            "x0": float(x[0]), "y1": float(y[-1]), "dx": dx, "nx": len(x), "ny": len(y),
            "bounds": [float(lon.min()), float(lat.min()), float(lon.max()), float(lat.max())], **stats}
    (out / f"{name}.json").write_text(json.dumps(meta))
    # Rebuild the index from every region on disk (safe when regions run in parallel); coarsest first.
    regions = [json.loads(p.read_text()) for p in out.glob("*.json") if p.name != "regions.json"]
    (out / "regions.json").write_text(json.dumps(sorted(regions, key=lambda r: -r["dx"]), indent=1))


def read(web: Path, name, lon, lat):
    """An exported region's ice surface (m, relative to today's ground) and ice mask at lon/lat points."""
    meta = json.loads((web / "data" / f"{name}.json").read_text())
    px = np.asarray(Image.open(web / "data" / f"{name}.png"))[:meta["ny"]].astype(np.float64)
    x, y = LAEA(meta["lon0"], meta["lat0"]).forward(lon, lat)
    rc = [(meta["y1"] - y) / meta["dx"], (x - meta["x0"]) / meta["dx"]]  # row 0 = north
    s = map_coordinates(px[..., 0] * 256 + px[..., 1] - 10000, rc, order=1, mode="nearest")
    return s, map_coordinates(px[..., 2], rc, order=0, mode="nearest") > 0
