"""Fit a smooth LGM temperature correction to the mapped ice: python -m glaciers.fit idaho [--iterations 8]

The climate models' LGM cooling has the wrong regional shape over the Rockies (too weak in the central
ranges, too strong near the ice sheets), so mapped glacier extents are used to correct it, as in climate
inversions from moraines. The correction lives on a 0.5-1° lon/lat grid of control points, so it can only shift
snowlines range by range, not draw individual glaciers. Held-out ranges are left out of the fit and test it.
Most iterations run on coarser grids (see model.levels); the last ones on the region's own grid.
"""

import argparse
import time

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates

from . import model, outlines
from .regions import correction_path, setup

HELD_OUT = {"idaho": {"Lemhi": (-113.60, -112.80, 44.00, 44.90), "Pioneer": (-114.25, -113.60, 43.50, 43.90)},
            "lost_river": {"Pioneer": (-114.25, -113.60, 43.50, 43.90)},
            "wind_river": {"East-side valleys": (-109.45, -108.95, 43.05, 43.50)},
            "sawtooth": {"White Clouds": (-114.75, -114.45, 44.00, 44.20)}}
SPACING = {"idaho": 1.0}  # degrees between control points; single ranges use 0.5, which also fits held-out ranges better
RADIUS = 50e3      # per degree of spacing: control points compare model and mapped ice within about this distance (m)
GAIN = 1.0         # K of correction per e-fold of area mismatch
MAX_STEP = 1.0     # K per iteration
LIMITS = (-4.0, 4.0)  # about the gap between climate models and glacier-based cooling estimates
SMOOTH = 0.5      # each iteration, blend half of the field with a smoothed copy (1 control point)
MIN_ICE = 100.0   # per degree² of spacing: km² of weighted ice before a control point responds
EDGE = 60e3       # ice within this distance of an inflow edge comes from outside and is not fitted (m)
SPINUP = 3000.0   # years of growth from nothing before the first iteration
ADJUST = 375.0    # years the ice adjusts per iteration on the region's own grid (longer on coarser ones)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("region", choices=HELD_OUT)
    p.add_argument("--iterations", type=int, default=8, help="at most this many per grid")
    a = p.parse_args()

    r = setup(a.region, use_correction=False)
    proj, spacing = r["proj"], SPACING.get(a.region, 0.5)
    radius, min_ice = RADIUS * spacing, MIN_ICE * spacing**2
    lons = np.arange(np.floor(r["lon"].min()), np.ceil(r["lon"].max()) + spacing, spacing)
    lats = np.arange(np.floor(r["lat"].min()), np.ceil(r["lat"].max()) + spacing, spacing)
    dT, H = np.zeros((len(lats), len(lons))), None
    for bed0, (T, P, z_ref), dx, edge in model.levels(r["bed0"], r["climate"], r["dx"], r["H_edge"]):
        lon, lat, X, Y = (model.resample(v, bed0.shape) for v in (r["lon"], r["lat"], *np.meshgrid(r["x"], r["y"])))
        mapped = model.resample(outlines.rasterize(proj, r["x"], r["y"]) * 1.0, bed0.shape) > 0.5
        box = lambda b: (lon > b[0]) & (lon < b[1]) & (lat > b[2]) & (lat < b[3])
        skip = np.zeros(bed0.shape, bool)  # cells left out of the fit: held-out ranges and edge inflow
        for b in HELD_OUT[a.region].values():
            skip |= box(b)
        if edge is not None:  # ice near an inflow edge reflects the parent run, not this climate
            skip |= distance_transform_edt(edge <= 0) * dx < EDGE
        weights = []  # per control point: (indices of nearby fitted cells, Gaussian weights)
        for px, py in (proj.forward(lo, la) for la in lats for lo in lons):
            d2 = (X - px) ** 2 + (Y - py) ** 2
            k = np.flatnonzero((d2 < (3 * radius) ** 2) & ~skip)
            weights.append((k, np.exp(-d2.ravel()[k] / (2 * radius**2))))
        overlap = lambda ice, where: (ice & mapped & where).sum() / max(((ice | mapped) & where).sum(), 1)

        H = None if H is None else model.refine(H, coarse, bed0)
        coarse = bed0
        for it in range(a.iterations):
            t0 = time.time()
            fitted = map_coordinates(dT, [(lat - lats[0]) / spacing, (lon - lons[0]) / spacing], order=1, mode="nearest")
            # From nothing the first time; after that the ice only has to adjust, faster on coarser grids.
            years = SPINUP if H is None else ADJUST * dx / r["dx"]
            H, _ = model.simulate(bed0, (T + fitted, P, z_ref), dx, years, edge, H)
            ice = H > 1

            # Update each control point from the weighted model and mapped ice areas around it.
            cell, step = dx * dx / 1e6, np.zeros_like(dT)
            for n, (k, w) in enumerate(weights):
                am, ao = (w * ice.ravel()[k]).sum() * cell, (w * mapped.ravel()[k]).sum() * cell
                if am + ao > min_ice:
                    step.flat[n] = np.clip(GAIN * np.log((am + min_ice / 2) / (ao + min_ice / 2)), -MAX_STEP, MAX_STEP)
            old, dT = dT, np.clip((1 - SMOOTH) * (dT + step) + SMOOTH * gaussian_filter(dT + step, 1.0, mode="nearest"), *LIMITS)
            np.savez(correction_path(a.region), lons=lons, lats=lats, dT=dT)

            held = "".join(f", held-out {name} {overlap(ice, box(b)):.0%}" for name, b in HELD_OUT[a.region].items())
            print(f"{dx:.0f} m grid, iteration {it + 1}: {time.time() - t0:.0f}s, fitted-area overlap {overlap(ice, ~skip):.0%}"
                  f"{held}, correction {dT.min():+.1f}..{dT.max():+.1f} K, changed by {np.abs(dT - old).max():.2f} K", flush=True)
            if np.abs(dT - old).max() < 0.1:
                break
    print(f"wrote {correction_path(a.region)}; now run: python -m glaciers {a.region}")


if __name__ == "__main__":
    main()
