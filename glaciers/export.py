"""Export a region's ice surface as one small PNG + metadata; the browser downscales it onto 30 m terrain."""

import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import distance_transform_edt, gaussian_filter


def write(web: Path, name, proj, x, y, dx, bed0, H, ocean, sea_level, stats):
    """PNG pixels: R,G = ice surface + 10000 (m, 16-bit), B = signed distance to the ice margin | 128 if exposed shelf.

    The surface is relative to today's ground (no isostatic depression), so the browser can compare it
    directly with present-day elevation tiles. The distance (grid cells, positive inside the ice, 1/4-cell steps,
    offset 64) is slightly smoothed so that its zero contour gives rounded margins rather than grid-cell steps,
    and lets the browser thin the ice towards the margin as real ice sheets do.
    """
    surface = np.clip(np.round(bed0 + H + 10000), 0, 65535).astype(np.uint32)
    shelf = ocean & (bed0 > sea_level)
    ice = H > 1
    # Margins are where ice ends on ground below its surface; ice ending against higher rock (nunataks,
    # valley walls) keeps its full thickness, so that rock doesn't count as outside.
    near = distance_transform_edt(~ice, return_distances=False, return_indices=True)
    outside = ~ice & (bed0 < (bed0 + H)[tuple(near)])
    sd = gaussian_filter(np.where(outside, 1 - distance_transform_edt(outside), distance_transform_edt(~outside)) - 0.5, 1.0)
    b = np.clip(np.round(sd * 4) + 64, 0, 127).astype(np.uint32) | (shelf * 128)
    rgb = np.stack([surface >> 8, surface & 255, b], -1).astype(np.uint8)
    out = web / "data"
    out.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgb[::-1]).save(out / f"{name}.png", optimize=True)  # row 0 = north

    lon, lat = proj.inverse(*np.meshgrid(x, y))
    meta = {"name": name, "lon0": float(np.degrees(proj.lon0)), "lat0": float(np.degrees(proj.lat0)),
            "x0": float(x[0]), "y1": float(y[-1]), "dx": dx, "nx": len(x), "ny": len(y),
            "bounds": [float(lon.min()), float(lat.min()), float(lon.max()), float(lat.max())], **stats}
    (out / f"{name}.json").write_text(json.dumps(meta))
    # Rebuild the index from every region on disk (safe when regions run in parallel); coarsest first.
    regions = [json.loads(p.read_text()) for p in out.glob("*.json") if p.name != "regions.json"]
    (out / "regions.json").write_text(json.dumps(sorted(regions, key=lambda r: -r["dx"]), indent=1))
