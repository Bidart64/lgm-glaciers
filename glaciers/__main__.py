"""python -m glaciers <region> [--dT -1] [--years N]"""

import argparse
import time

import numpy as np

from . import export, model, outlines
from .regions import REGIONS, WEB, setup


def main():
    p = argparse.ArgumentParser(description="Compute Last Glacial Maximum glaciers.")
    p.add_argument("region", choices=REGIONS)
    p.add_argument("--dT", type=float, help="extra temperature offset on top of PMIP4 LGM climate in K (default per region)")
    p.add_argument("--sea-level", type=float, default=-120.0, help="LGM sea level in m")
    p.add_argument("--years", type=float, help="simulated years (default per region)")
    a = p.parse_args()

    r = setup(a.region, dT=a.dT)
    dx, x, y = r["dx"], r["x"], r["y"]
    t0 = time.time()
    print(f"{a.region}: {len(y)}x{len(x)} grid, {dx / 1000:g} km, extra dT={r['dT']} K"
          f"{', fitted correction' if r['corrected'] else ''}")
    H, bed = model.simulate(r["bed0"], r["climate"], dx, a.years or r["years"], r["dt"], a.sea_level, r["H_edge"])
    s = np.maximum(bed + H, a.sea_level)
    flow = (*model.surface_velocity(H, s, dx, model.sliding(r["bed0"], dx)), model.mass_balance(s, *r["climate"][:3]))
    area = (H > 1).sum() * dx * dx / 1e6
    mapped = outlines.rasterize(r["proj"], x, y)
    iou = ((H > 1) & mapped).sum() / max(((H > 1) | mapped).sum(), 1)  # agreement with field evidence
    print(f"done in {time.time() - t0:.0f}s: ice {area:,.0f} km2, max {H.max():,.0f} m thick")
    print(f"mapped LGM ice {mapped.sum() * dx * dx / 1e6:,.0f} km2, overlap (intersection / union) {iou:.2f}")

    export.write(WEB, a.region, r["proj"], x, y, dx, r["bed0"], H, r["climate"][3], a.sea_level,
                 {"dT": r["dT"], "corrected": r["corrected"], "sea_level": a.sea_level, "ice_km2": round(area),
                  "max_ice_m": round(float(H.max())), "mapped_km2": round(mapped.sum() * dx * dx / 1e6),
                  "iou": round(float(iou), 2)}, flow)
    print(f"wrote {WEB / 'data' / a.region}.png\nview:  python -m http.server -d web  then open http://localhost:8000")


if __name__ == "__main__":
    main()
