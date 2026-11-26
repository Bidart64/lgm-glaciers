"""Fit a smooth LGM temperature correction to the mapped ice: python -m glaciers.fit sawtooth_range [--iterations 8] [--resume]

The climate models' LGM cooling has the wrong regional shape over the Rockies (too weak in the central
ranges, too strong near the ice sheets), so mapped glacier extents are used to correct it, as in climate
inversions from moraines. The correction lives on a 0.5-1° lon/lat grid of control points, so it can only shift
snowlines range by range, not draw individual glaciers. Held-out ranges are left out of the fit and test it.
Most iterations run on coarser grids (see model.levels); the last three on the region's own grid.
"""

import argparse
import time

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter, map_coordinates

from . import model
from .regions import PATCH_NAMES, REGIONS, correction_path, inflow_scale, setup

HELD_OUT = {"alaska": {"Central Brooks Range": (-154.0, -148.0, 67.4, 68.8)}}  # name: (west, east, south, north)


def held_out(region, lon, lat):
    """{name: cells} of the ground left out of a region's fit to test it: named boxes where there are any, or
    else one in five quarter-degree squares, scattered."""
    if region in HELD_OUT:
        return {name: (lon > w) & (lon < e) & (lat > s) & (lat < n) for name, (w, e, s, n) in HELD_OUT[region].items()}
    return {"squares": (np.floor(lon * 4) + 2 * np.floor(lat * 4)) % 5 == 0}


SPACING = {"alaska": 1.0}  # degrees between control points; single ranges use 0.5, which also fits held-out ranges better
RADIUS = 50e3      # per degree of spacing: control points compare model and mapped ice within about this distance (m)
GAIN = 1.0         # K of correction per e-fold of area mismatch
MAX_STEP = 1.0     # K per iteration
LIMITS = (-4.0, 4.0)  # about the gap between climate models and glacier-based cooling estimates
SMOOTH = 0.5      # each iteration, blend half of the field with a smoothed copy (1 control point)
MIN_ICE = 100.0   # per degree² of spacing: km² of weighted ice before a control point responds
EDGE = 60e3       # ice within this distance of an inflow edge comes from outside and is not fitted (m)
INFLOW_GAIN = 0.5          # how much of the lobe's excess (as a log ratio of areas) each iteration takes off the inflow
INFLOW_LIMITS = 0.1, 1.5   # shares of the coarse run's edge thickness
FINAL = 3         # iterations on the region's own grid, where each one costs a full run
TEST_AREA = 500.0  # km² of mapped ice below which nothing is held out: a test square could swallow the whole range


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("region", choices=[r for r in REGIONS if r != "north_america"])
    p.add_argument("--iterations", type=int, default=8, help="at most this many per coarser grid")
    p.add_argument("--resume", action="store_true", help="start from the region's existing correction")
    a = p.parse_args()

    r = setup(a.region, use_correction=False)  # the climate as the models give it; the fit adds dT to it
    proj, spacing = r["proj"], SPACING.get(a.region, 0.5)
    area = r["mapped"].sum() * r["dx"] ** 2 / 1e6  # km² of mapped ice
    radius, min_ice = RADIUS * spacing, min(MIN_ICE * spacing**2, area / 10)  # small ranges must still respond
    lons = np.arange(np.floor(r["lon"].min()), np.ceil(r["lon"].max()) + spacing, spacing)
    lats = np.arange(np.floor(r["lat"].min()), np.ceil(r["lat"].max()) + spacing, spacing)
    dT = np.zeros((len(lats), len(lons)))  # the correction: one value per control point (K)
    inflow = 1.0  # and, for a patch the ice sheet flows into, the share of the coarse run's thickness let in
    if a.resume and correction_path(a.region).exists():  # carry on from an earlier fit, which converges sooner
        dT, inflow = np.load(correction_path(a.region))["dT"], inflow_scale(a.region)
    for bed0, (T, P, z_ref), dx, edge in model.levels(r["bed0"], r["climate"], r["dx"], r["H_edge"]):
        lon, lat, X, Y = (model.resample(v, bed0.shape) for v in (r["lon"], r["lat"], *np.meshgrid(r["x"], r["y"])))
        mapped = model.resample(r["mapped"] * 1.0, bed0.shape) > 0.5
        tests = held_out(a.region, lon, lat) if area > TEST_AREA else {}
        # Cells left out of the fit: where the outlines say nothing, held-out ranges and edge inflow.
        known = model.resample(r["known"] * 1.0, bed0.shape) > 0.5
        skip = ~known
        for cells in tests.values():
            skip |= cells
        lobe = np.zeros_like(skip)  # ground near an inflow edge, whose ice reflects the inflow, not this climate
        if edge is not None:
            lobe = (distance_transform_edt(edge <= 0) * dx < EDGE) & ~skip
            skip |= lobe
        weights = []  # per control point: (indices of nearby fitted cells, Gaussian weights)
        for px, py in (proj.forward(lo, la) for la in lats for lo in lons):
            d2 = (X - px) ** 2 + (Y - py) ** 2
            k = np.flatnonzero((d2 < (3 * radius) ** 2) & ~skip)
            weights.append((k, np.exp(-d2.ravel()[k] / (2 * radius**2))))
        overlap = lambda ice, where: (ice & mapped & where).sum() / max(((ice | mapped) & where).sum(), 1)

        gain, last = np.ones_like(dT), np.zeros_like(dT)  # per control point: how much of its step is taken
        for it in range(a.iterations if dx > r["dx"] else FINAL):
            t0 = time.time()
            fitted = map_coordinates(dT, [(lat - lats[0]) / spacing, (lon - lons[0]) / spacing], order=1, mode="nearest")
            # Grown from nothing every time, exactly as a run does: with ice left from colder iterations the
            # fit would describe a state that a run never reaches.
            H, _ = model.grow(bed0, (T + fitted, P, z_ref), dx, r["years"], None if edge is None else edge * inflow)
            ice = H > model.GLACIER

            # Update each control point from the weighted model and mapped ice areas around it.
            cell, step = dx * dx / 1e6, np.zeros_like(dT)
            for n, (k, w) in enumerate(weights):
                am, ao = (w * ice.ravel()[k]).sum() * cell, (w * mapped.ravel()[k]).sum() * cell  # model, mapped (km²)
                if am + ao > min_ice:  # too much model ice here: warm; too little: cool
                    step.flat[n] = np.clip(GAIN * np.log((am + min_ice / 2) / (ao + min_ice / 2)), -MAX_STEP, MAX_STEP)
            # Where a basin is close to filling with ice, a small change flips it between valley glaciers and an
            # ice cap, and the fit swings back and forth: a point whose step reverses takes half as much from then on.
            gain[step * last < 0] /= 2
            last, step = step, step * gain
            # Take the step, then pull each point part of the way towards its neighbours and keep it in bounds.
            old, dT = dT, np.clip((1 - SMOOTH) * (dT + step) + SMOOTH * gaussian_filter(dT + step, 1.0, mode="nearest"), *LIMITS)
            # The coarse run's ice sheet is too thick at the patches' edges and its lobes run too far: the inflow
            # is scaled until the lobe covers as much ground as the mapped one.
            was = inflow
            if a.region in PATCH_NAMES and lobe.any():
                excess = np.log(((ice & lobe).sum() + 1) / ((mapped & lobe).sum() + 1))
                inflow = float(np.clip(inflow * np.exp(-INFLOW_GAIN * np.clip(excess, -1, 1)), *INFLOW_LIMITS))
            score = overlap(ice, ~skip)
            np.savez(correction_path(a.region), lons=lons, lats=lats, dT=dT, inflow=inflow)

            held = "".join(f", held-out {name} {overlap(ice, cells & known):.0%}" for name, cells in tests.items())
            print(f"{dx:.0f} m grid, iteration {it + 1}: {time.time() - t0:.0f}s, fitted-area overlap {score:.0%}"
                  f"{held}, correction {dT.min():+.1f}..{dT.max():+.1f} K, changed by {np.abs(dT - old).max():.2f} K"
                  + (f", lobe overlap {overlap(ice, lobe):.0%}, inflow x{inflow:.2f}" if lobe.any() else ""), flush=True)
            if np.abs(dT - old).max() < 0.1 and abs(np.log(inflow / was)) < 0.05:
                break
    print(f"wrote {correction_path(a.region)}; now run: python -m glaciers {a.region}")


if __name__ == "__main__":
    main()
