"""python -m glaciers <region> [--years N] [--mapped]"""

import argparse
import time

import numpy as np

from . import export, model
from .regions import REGIONS, WEB, setup


def main():
    p = argparse.ArgumentParser(description="Compute Last Glacial Maximum glaciers.")
    p.add_argument("region", choices=REGIONS)
    p.add_argument("--years", type=float, help="simulated years (default per region)")
    p.add_argument("--mapped", action="store_true", help="fill the mapped extent with ice instead of predicting it")
    a = p.parse_args()

    r = setup(a.region, mapped=a.mapped)
    dx, x, y, bed0 = r["dx"], r["x"], r["y"], r["bed0"]
    t0 = time.time()
    print(f"{a.region}: {len(y)}x{len(x)} grid, {dx / 1000:g} km{', fitted correction' if r['corrected'] else ''}")
    H, bed = model.grow(bed0, r["climate"], dx, a.years or r["years"], r["H_edge"])
    s = np.maximum(bed + H, model.SEA_LEVEL)
    vx, vy = model.surface_velocity(H, s, dx, model.sliding(bed0, dx))
    ice, mapped = H > model.GLACIER, r["mapped"]
    km2 = lambda mask: round(mask.sum() * dx * dx / 1e6)
    known = r["known"]  # agreement with field evidence, where there is any
    iou = (ice & mapped & known).sum() / max(((ice | mapped) & known).sum(), 1)
    print(f"done in {time.time() - t0:.0f}s: ice {km2(ice):,} km2, max {H.max():,.0f} m thick")
    print(f"mapped LGM ice {km2(mapped):,} km2, overlap (intersection / union) {iou:.2f}")

    export.write(WEB, a.region, r["proj"], x, y, dx, bed0, r["rebuilt"], H, vx, vy, model.mass_balance(s, *r["natural"]),
                 {"corrected": r["corrected"], "mapped": a.mapped, "ice_km2": km2(ice), "max_ice_m": round(float(H.max())),
                  "mapped_km2": km2(mapped), "iou": round(float(iou), 2)})
    print(f"wrote {WEB / 'data' / a.region}.png\nview:  python -m http.server -d web  then open http://localhost:8000")


if __name__ == "__main__":
    main()
