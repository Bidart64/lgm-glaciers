"""python -m glaciers <region> [--dT -1] [--years N]"""

import argparse
import time
from pathlib import Path

import numpy as np

from . import dem, export, model, outlines, pmip

# name: (centre lon, centre lat, half-width km, half-height km, grid km, years, time step yr, extra cooling K)
# The 2-3° PMIP4 models resolve the ice sheets but not mountain ranges, whose LGM cooling they underestimate.
REGIONS = {
    "alps":          (9.5, 46.6, 600, 380, 2, 3000, 10, -3.5),
    "pyrenees":      (0.5, 42.6, 300, 150, 2, 3000, 10, -3.5),
    "new_zealand":   (171.0, -43.0, 400, 450, 3, 3000, 10, -3.5),
    "patagonia":     (-71.0, -47.0, 450, 1100, 4, 5000, 10, -3.5),
    "europe":        (20.0, 65.0, 2800, 2200, 20, 40000, 25, 0),
    "north_america": (-95.0, 58.0, 3300, 2600, 20, 40000, 25, 0),
}


def main():
    p = argparse.ArgumentParser(description="Compute Last Glacial Maximum glaciers.")
    p.add_argument("region", choices=REGIONS)
    p.add_argument("--dT", type=float, help="extra temperature offset on top of PMIP4 LGM climate in K (default per region)")
    p.add_argument("--sea-level", type=float, default=-120.0, help="LGM sea level in m")
    p.add_argument("--years", type=float, help="simulated years (default per region)")
    a = p.parse_args()

    lon0, lat0, hw, hh, km, years, dt, dT = REGIONS[a.region]
    dT = dT if a.dT is None else a.dT
    dx = km * 1000.0
    proj = dem.LAEA(lon0, lat0)
    x = np.arange(-hw * 1000, hw * 1000 + 1, dx)
    y = np.arange(-hh * 1000, hh * 1000 + 1, dx)
    bed0 = dem.load_bed(proj, x, y, dx)
    T, P, z_ref, ocean = dem.load_climate(proj, x, y)
    lgm_dT, ratio = pmip.sample(*proj.inverse(*np.meshgrid(x, y)))
    climate = T + lgm_dT + dT, P * ratio, z_ref, ocean  # present climate + LGM anomalies

    t0 = time.time()
    print(f"{a.region}: {len(y)}x{len(x)} grid, {dx / 1000:g} km, extra dT={dT} K")
    H, _ = model.simulate(bed0, climate, dx, a.years or years, dt, a.sea_level)
    area = (H > 1).sum() * dx * dx / 1e6
    mapped = outlines.rasterize(proj, x, y)
    iou = ((H > 1) & mapped).sum() / max(((H > 1) | mapped).sum(), 1)  # agreement with field evidence
    print(f"done in {time.time() - t0:.0f}s: ice {area:,.0f} km2, max {H.max():,.0f} m thick")
    print(f"mapped LGM ice {mapped.sum() * dx * dx / 1e6:,.0f} km2, overlap (intersection / union) {iou:.2f}")

    web = Path(__file__).resolve().parent.parent / "web"
    export.write(web, a.region, proj, x, y, dx, bed0, H, ocean, a.sea_level,
                 {"dT": dT, "sea_level": a.sea_level, "ice_km2": round(area), "max_ice_m": round(float(H.max())),
                  "mapped_km2": round(mapped.sum() * dx * dx / 1e6), "iou": round(float(iou), 2)})
    print(f"wrote {web / 'data' / a.region}.png\nview:  python -m http.server -d web  then open http://localhost:8000")


if __name__ == "__main__":
    main()
