"""Fit a smooth LGM temperature correction to the mapped ice: python -m glaciers.fit idaho [--iterations 8]

The climate models' LGM cooling has the wrong regional shape over Idaho (too weak in the central ranges,
too strong near the ice sheets), so mapped glacier extents are used to correct it, as in climate inversions
from moraines. The correction lives on a 1° lon/lat grid of control points, so it can only shift snowlines
region by region, not draw individual glaciers. Held-out ranges are left out of the fit and test it.
"""

import argparse
import time

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter

from . import model, outlines
from .regions import correction, correction_path, setup

HELD_OUT = {"idaho": {"Lemhi": (-113.60, -112.80, 44.00, 44.90), "Pioneer": (-114.25, -113.60, 43.50, 43.90)},
            "wind_river": {"East-side valleys": (-109.45, -108.95, 43.05, 43.50)},
            "sawtooth": {"White Clouds": (-114.75, -114.45, 44.00, 44.20)}}
RADIUS = 50e3      # control points compare model and mapped ice within about this distance (m)
GAIN = 1.0         # K of correction per e-fold of area mismatch
MAX_STEP = 1.0     # K per iteration
LIMITS = (-4.0, 4.0)  # about the gap between climate models and glacier-based cooling estimates
SMOOTH = 0.5      # each iteration, blend half of the field with a smoothed copy (1 control point)
MIN_ICE = 100.0   # km² of weighted ice before a control point responds; smaller patches are below the grid
EDGE = 60e3       # ice within this distance of an inflow edge comes from outside and is not fitted (m)
FIT_KM = {"idaho": 1.0, "wind_river": 0.5, "sawtooth": 0.2}  # fitting grid (km); the correction itself is resolution-independent


def save(region, lons, lats, dT, shift):
    """The correction, plus the PRISM baseline shift it was fitted with (reused by regions that borrow it)."""
    extra = {} if shift is None else {"shift": shift}
    np.savez(correction_path(region), lons=lons, lats=lats, dT=dT, **extra)


def box(lon, lat, b):
    return (lon > b[0]) & (lon < b[1]) & (lat > b[2]) & (lat < b[3])


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("region", choices=HELD_OUT)
    p.add_argument("--iterations", type=int, default=8)
    a = p.parse_args()

    r = setup(a.region, km=FIT_KM[a.region], use_correction=False)
    lon, lat, dx, proj = r["lon"], r["lat"], r["dx"], r["proj"]
    mapped = outlines.rasterize(proj, r["x"], r["y"])
    skip = np.zeros(lon.shape, bool)  # cells left out of the fit: held-out ranges and edge inflow
    if r["H_edge"] is not None:  # ice near an inflow edge reflects the parent run, not this climate
        inflow = np.zeros(lon.shape, bool)
        inflow[r["H_edge"] > 0] = True
        skip |= distance_transform_edt(~inflow) * dx < EDGE
    for b in HELD_OUT[a.region].values():
        skip |= box(lon, lat, b)

    lons = np.arange(np.floor(lon.min()), np.ceil(lon.max()) + 1)
    lats = np.arange(np.floor(lat.min()), np.ceil(lat.max()) + 1)
    dT = np.zeros((len(lats), len(lons)))
    X, Y = np.meshgrid(r["x"], r["y"])
    weights = []  # per control point: (indices of nearby fitted cells, Gaussian weights)
    for la in lats:
        for lo in lons:
            px, py = proj.forward(lo, la)
            d2 = (X - px) ** 2 + (Y - py) ** 2
            k = np.flatnonzero((d2 < (3 * RADIUS) ** 2) & ~skip)
            weights.append((k, np.exp(-d2.ravel()[k] / (2 * RADIUS**2))))

    T0, P, z_ref, ocean = r["climate"]
    H = None
    for it in range(a.iterations):
        save(a.region, lons, lats, dT, r["prism_shift"])  # correction() reads it back
        climate = T0 + correction(a.region, lon, lat), P, z_ref, ocean
        t0 = time.time()
        H, _ = model.simulate(r["bed0"], climate, dx, 3000 if H is None else 1500, r["dt"], -120.0, r["H_edge"], H)
        ice = H > 1

        # Update each control point from the weighted model and mapped ice areas around it.
        cell = dx * dx / 1e6
        step = np.zeros_like(dT)
        for n, (k, w) in enumerate(weights):
            am, ao = (w * ice.ravel()[k]).sum() * cell, (w * mapped.ravel()[k]).sum() * cell
            if am + ao > MIN_ICE:
                step.flat[n] = np.clip(GAIN * np.log((am + MIN_ICE / 2) / (ao + MIN_ICE / 2)), -MAX_STEP, MAX_STEP)
        dT = dT + step
        dT = np.clip((1 - SMOOTH) * dT + SMOOTH * gaussian_filter(dT, 1.0, mode="nearest"), *LIMITS)

        fit = (ice & mapped & ~skip).sum() / max(((ice | mapped) & ~skip).sum(), 1)
        line = [f"iteration {it + 1}: {time.time() - t0:.0f}s, fitted-area overlap {fit:.0%}"]
        for name, b in HELD_OUT[a.region].items():
            k = box(lon, lat, b)
            line.append(f"held-out {name} {(ice & mapped & k).sum() / max(((ice | mapped) & k).sum(), 1):.0%}")
        line.append(f"correction {dT.min():+.1f}..{dT.max():+.1f} K, largest step {np.abs(step).max():.2f} K")
        print(", ".join(line), flush=True)
        if np.abs(step).max() < 0.2:
            break
    save(a.region, lons, lats, dT, r["prism_shift"])
    print(f"wrote {correction_path(a.region)}; now run: python -m glaciers {a.region}")


if __name__ == "__main__":
    main()
