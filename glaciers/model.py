"""Last Glacial Maximum ice: degree-day mass balance + shallow-ice flow + bed sinking."""

import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter, map_coordinates, uniform_filter
from scipy.sparse.linalg import spsolve
from scipy.special import erfc

RHO_I, RHO_M = 917.0, 3300.0
G = 9.81
A = 1e-16          # Glen flow-law rate factor (Pa^-3 yr^-1)
N = 3
GAMMA = 2 * A * (RHO_I * G) ** N / (N + 2)
LAPSE = 0.0065     # K per m
SIGMA = 3.0        # daily temperature variability around the monthly mean (K)
DDF_SNOW, DDF_ICE = 0.003, 0.008  # melt per positive degree-day (m water / K / day)
REFREEZE = 0.6     # fraction of snowmelt that refreezes in the snowpack
SEA_LEVEL = -120.0 # LGM sea level (m)
MAX_WATER = 500.0  # ice calves rather than grounds in deeper water (m)
TAU_ISO = 3000.0   # bed relaxation time (yr)
SLIDE = 3e-3       # basal sliding: velocity per basal shear stress (m / yr / Pa)
SOFT_BED = 5.0     # sliding enhancement on smooth lowland (sediment) beds
STEP_SPEED = 80.0  # time step = grid spacing / this (m/yr): results stop changing below it
MAX_STEP = 25.0    # yr
ADJUST = 1 / 6     # share of a run's years spent adjusting on the full grid (twice that on the 2x coarser one)


def mass_balance(s, T_ref, P, z_ref):
    """Annual surface mass balance (m ice / yr) from monthly climate, lapse-corrected to surface s."""
    T = T_ref + LAPSE * (z_ref - s)
    snow = (P * np.clip((2 - T) / 4, 0, 1)).sum(0) / 1000  # rain/snow mix between -2 and 2 °C (m water)
    pdd = (30.4 * (SIGMA / np.sqrt(2 * np.pi) * np.exp(-T**2 / (2 * SIGMA**2))
                   + T / 2 * erfc(-T / (np.sqrt(2) * SIGMA)))).sum(0)
    melt = (1 - REFREEZE) * np.minimum(pdd * DDF_SNOW, snow) + np.maximum(pdd - snow / DDF_SNOW, 0) * DDF_ICE
    return (snow - melt) * 1000 / RHO_I


def sliding(bed0, dx):
    """Basal sliding coefficient per cell: ice slides fast over smooth, low-lying sediment (e.g. the
    Laurentide's southern lobes)."""
    mean = gaussian_filter(bed0, 50e3 / dx)
    rough = np.sqrt(np.maximum(gaussian_filter(bed0**2, 50e3 / dx) - mean**2, 0))
    return SLIDE * np.where((rough < 150) & (bed0 < 800), SOFT_BED, 1.0)


def surface_velocity(H, s, dx, slide):
    """Shallow-ice surface velocity (m/yr, x east and y north): deformation plus sliding, downhill."""
    gy, gx = np.gradient(s, dx)
    slope = np.maximum(np.hypot(gx, gy), 1e-9)
    speed = 2 * A / (N + 1) * (RHO_I * G * slope) ** N * H ** (N + 1) + slide * RHO_I * G * H * slope
    return -speed * gx / slope, -speed * gy / slope


def _faces(s, H, slide, dx, dt):
    """Along axis 1, on the faces between cells: the ice flux (m²/yr, towards higher index) and what a
    Newton step needs: its sensitivity to the surface slope and to the thickness of the upstream cell
    (split by whether that is the lower- or higher-index cell), both scaled for the matrix."""
    grad = np.diff(s, axis=1) / dx
    up = grad > 0  # ice flows down from the higher-index cell
    # Upstream thickness: an ice-free cell above a glacier (a valley wall) gives no ice. Averaging the two
    # cells instead makes such cells go negative, and clipping that creates ice.
    Hu = np.where(up, H[:, 1:], H[:, :-1])
    deform = GAMMA * Hu ** (N + 2) * grad**2
    slip = 0.5 * (slide[:, 1:] + slide[:, :-1]) * RHO_I * G * Hu**2
    a = dt / dx * np.abs(((N + 2) * deform + 2 * slip) * grad) / np.maximum(Hu, 1e-9)
    return -(deform + slip) * grad, dt / dx**2 * (N * deform + slip), np.where(up, 0, a), np.where(up, a, 0)


def simulate(bed0, climate, dx, years, H_edge=None, H0=None):
    """Grow ice on `bed0` under LGM `climate` for `years` and return (ice thickness, depressed bed).

    `H_edge` holds the domain edges at a parent run's ice thickness, so ice can flow in from outside.
    `H0` starts from existing ice instead of none (with the bed already sunk under it).
    """
    slide = sliding(bed0, dx)
    H = np.zeros_like(bed0) if H0 is None else H0.copy()
    edge = np.zeros(H.shape, bool)
    edge[[0, -1]] = edge[:, [0, -1]] = True
    bed = bed0 - RHO_I / RHO_M * H
    ny, nx = H.shape
    dt = min(dx / STEP_SPEED, MAX_STEP)
    for _ in range(int(years / dt)):
        s = np.maximum(bed + H, SEA_LEVEL)
        M = mass_balance(s, *climate)
        melt = np.maximum(M, -H / dt)  # in the solve, melt only the ice that is there; the rest comes after
        # One Newton step of the implicit equation H' = H + dt (M - div flux(H')). The flux depends on H' both
        # through the slope (diffusion) and through the upstream thickness (a wave moving with the ice);
        # holding the second at its old value, as simpler schemes do, is unstable at these time steps and
        # the oscillations create ice.
        fx, dxx, ax_lo, ax_hi = _faces(s, H, slide, dx, dt)
        fy, dyy, ay_lo, ay_hi = (v.T for v in _faces(s.T, H.T, slide.T, dx, dt))
        diag, rhs = np.ones((ny, nx)), dt * melt
        east, west, north, south = (np.zeros((ny, nx)) for _ in range(4))
        diag[:, :-1] += dxx + ax_lo; diag[:, 1:] += dxx + ax_hi; diag[:-1] += dyy + ay_lo; diag[1:] += dyy + ay_hi
        east[:, :-1], west[:, 1:], north[:-1], south[1:] = -dxx - ax_hi, -dxx - ax_lo, -dyy - ay_hi, -dyy - ay_lo
        rhs[:, :-1] -= dt / dx * fx; rhs[:, 1:] += dt / dx * fx; rhs[:-1] -= dt / dx * fy; rhs[1:] += dt / dx * fy
        mat = sparse.diags([diag.ravel(), east.ravel()[:-1], west.ravel()[1:], north.ravel()[:-nx], south.ravel()[nx:]],
                           [0, 1, -1, nx, -nx], format="csc")
        # A direct sparse solve: faster here than multigrid-preconditioned iterations, which converge
        # poorly because fast ice-sheet flow and slow mountain ice make the matrix ill-conditioned.
        H = np.maximum(H + spsolve(mat, rhs.ravel()).reshape(ny, nx) + dt * (M - melt), 0.0)
        H[bed < SEA_LEVEL - MAX_WATER] = 0.0  # ice calves away over deep water
        H[edge] = 0.0 if H_edge is None else H_edge[edge]
        bed += dt / TAU_ISO * (bed0 - RHO_I / RHO_M * H - bed)
    return H, bed


def resample(a, shape):
    """Bilinear resampling of the last two axes to `shape` (box-averaged first when coarsening)."""
    if a.ndim > 2:
        return np.array([resample(m, shape) for m in a])
    if shape[0] < a.shape[0]:
        a = uniform_filter(a, round(a.shape[0] / shape[0]), mode="nearest")
    rc = np.meshgrid(*(np.linspace(0, n - 1, m) for n, m in zip(a.shape, shape)), indexing="ij")
    return map_coordinates(a, rc, order=1)


def levels(bed0, climate, dx, H_edge=None):
    """The same problem on grids 4x and 2x coarser and on the full grid: (bed0, climate, dx, H_edge) each.
    Coarse grids settle the ice cheaply; each finer one starts from the last and only has to adjust.
    Grids already coarser than 2 km are used as they are."""
    for f in (4, 2, 1) if dx < 2000 else (1,):
        shape = tuple((n - 1) // f + 1 for n in bed0.shape)
        c = lambda a: a if f == 1 or a is None else resample(a, shape)
        yield c(bed0), tuple(c(a) for a in climate), dx * (bed0.shape[0] - 1) / (shape[0] - 1), c(H_edge)


def refine(H, bed0, fine_bed0):
    """Ice from a coarser grid on a finer one: its surface, not its thickness, is what carries over."""
    s = resample(bed0 + H, fine_bed0.shape)
    return np.where(resample((H > 0) * 1.0, fine_bed0.shape) > 0.5, np.maximum(s - fine_bed0, 0), 0)


def grow(bed0, climate, dx, years, H_edge=None):
    """Ice grown from nothing for `years` on the coarsest grid, then adjusted on each finer one (see levels).
    Returns (ice thickness, depressed bed)."""
    H = None
    for bed0_l, climate_l, dx_l, edge in levels(bed0, climate, dx, H_edge):
        H0 = None if H is None else refine(H, coarse, bed0_l)
        H, bed = simulate(bed0_l, climate_l, dx_l, years if H is None else years * ADJUST * dx_l / dx, edge, H0)
        coarse = bed0_l
    return H, bed
