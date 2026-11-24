"""Export a region's ice surface as one small PNG + metadata; the browser downscales it onto 30 m terrain."""

import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates

from .dem import LAEA


def write(web: Path, name, proj, x, y, dx, bed0, H, ocean, sea_level, stats, flow=None):
    """PNG pixels: R,G = ice surface + 10000 (m, 16-bit; next to the ice, the nearest ice surface),
    B = signed distance to the ice margin | 128 if exposed shelf.

    The surface is relative to today's ground (no isostatic depression), so the browser can compare it
    directly with present-day elevation tiles. The distance (grid cells, positive inside the ice, 1/4-cell steps,
    offset 64) is slightly smoothed so that its zero contour gives rounded margins rather than grid-cell steps,
    and lets the browser thin the ice towards the margin as real ice sheets do.
    """
    shelf = ocean & (bed0 > sea_level)
    ice = H > 1
    dist, near = distance_transform_edt(~ice, return_indices=True)
    ice_surface = (bed0 + H)[tuple(near)]  # surface of the nearest ice
    # Just outside the ice (valley walls), carry the ice surface rather than the rock, so the browser's
    # interpolated surface stays flat across a valley and ends where it meets the walls, instead of
    # ramping up them.
    surface = np.where(ice | (dist > 2), bed0 + H, ice_surface)
    surface = np.clip(np.round(surface + 10000), 0, 65535).astype(np.uint32)
    # Margins are where ice ends on ground below its surface; ice ending against higher rock (nunataks,
    # valley walls) keeps its full thickness, so that rock doesn't count as outside. Beyond two cells from
    # the ice everything is outside.
    outside = ~ice & ((bed0 < ice_surface) | (dist > 2))
    sd = gaussian_filter(np.where(outside, 1 - distance_transform_edt(outside), distance_transform_edt(~outside)) - 0.5, 1.0)
    b = np.clip(np.round(sd * 4) + 64, 0, 127).astype(np.uint32) | (shelf * 128)
    rgb = np.stack([surface >> 8, surface & 255, b], -1).astype(np.uint8)
    out = web / "data"
    out.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgb[::-1]).save(out / f"{name}.png", optimize=True)  # row 0 = north
    if flow is not None:
        Image.fromarray(surface_detail(dx, H, *flow)[::-1]).save(out / f"{name}.flow.png", optimize=True)

    lon, lat = proj.inverse(*np.meshgrid(x, y))
    meta = {"name": name, "lon0": float(np.degrees(proj.lon0)), "lat0": float(np.degrees(proj.lat0)),
            "x0": float(x[0]), "y1": float(y[-1]), "dx": dx, "nx": len(x), "ny": len(y),
            "bounds": [float(lon.min()), float(lat.min()), float(lon.max()), float(lat.max())], **stats}
    (out / f"{name}.json").write_text(json.dumps(meta))
    # Rebuild the index from every region on disk (safe when regions run in parallel); coarsest first.
    regions = [json.loads(p.read_text()) for p in out.glob("*.json") if p.name != "regions.json"]
    (out / "regions.json").write_text(json.dumps(sorted(regions, key=lambda r: -r["dx"]), indent=1))


CREVASSE_STRAIN = 0.02  # 1/yr: tensile strain rate at which crevasses are fully developed


def surface_detail(dx, H, vx, vy, M):
    """What the ice surface looks like, as a second PNG for the viewer:
    R = crevassing (0-255, from the largest tensile strain rate), G = direction of that stretching
    (0-255 for 0-180°, east = 0, anticlockwise; crevasses open across it), B = 128 + 25 × mass balance (m/yr),
    so above 128 is snow-covered accumulation area and below is bare ice."""
    dvx_dy, dvx_dx = np.gradient(vx, dx)
    dvy_dy, dvy_dx = np.gradient(vy, dx)
    exy = 0.5 * (dvx_dy + dvy_dx)
    half = 0.5 * (dvx_dx - dvy_dy)
    e1 = 0.5 * (dvx_dx + dvy_dy) + np.hypot(half, exy)  # largest principal strain rate
    angle = 0.5 * np.arctan2(exy, half) % np.pi           # its direction
    crevasse = np.where(H > 1, np.clip(e1 / CREVASSE_STRAIN, 0, 1), 0)
    return np.stack([np.round(crevasse * 255), np.round(angle / np.pi * 255) % 256,
                     np.clip(np.round(128 + 25 * M), 0, 255)], -1).astype(np.uint8)


def read(web: Path, name, lon, lat):
    """An exported region's ice surface (m, relative to today's ground) and ice mask at lon/lat points."""
    meta = json.loads((web / "data" / f"{name}.json").read_text())
    px = np.asarray(Image.open(web / "data" / f"{name}.png")).astype(np.float64)
    x, y = LAEA(meta["lon0"], meta["lat0"]).forward(lon, lat)
    rc = [(meta["y1"] - y) / meta["dx"], (x - meta["x0"]) / meta["dx"]]  # row 0 = north
    s = map_coordinates(px[..., 0] * 256 + px[..., 1] - 10000, rc, order=1, mode="nearest")
    d = map_coordinates((px[..., 2] % 128) - 64, rc, order=1, mode="nearest")
    return s, d > 0
