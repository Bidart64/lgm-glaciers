"""Last Glacial Maximum ice: degree-day mass balance + shallow-ice flow + bed sinking."""

import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter
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
MAX_WATER = 500.0  # ice calves rather than grounds in deeper water (m)
TAU_ISO = 3000.0   # bed relaxation time (yr)
SLIDE = 3e-3       # basal sliding: velocity per basal shear stress (m / yr / Pa)
SOFT_BED = 5.0     # sliding enhancement on smooth lowland (sediment) beds


def mass_balance(s, T_ref, P, z_ref):
    """Annual surface mass balance (m ice / yr) from monthly climate, lapse-corrected to surface s."""
    T = T_ref + LAPSE * (z_ref - s)
    snow = (P * np.clip((2 - T) / 4, 0, 1)).sum(0) / 1000  # rain/snow mix between -2 and 2 °C (m water)
    pdd = (30.4 * (SIGMA / np.sqrt(2 * np.pi) * np.exp(-T**2 / (2 * SIGMA**2))
                   + T / 2 * erfc(-T / (np.sqrt(2) * SIGMA)))).sum(0)
    melt = (1 - REFREEZE) * np.minimum(pdd * DDF_SNOW, snow) + np.maximum(pdd - snow / DDF_SNOW, 0) * DDF_ICE
    return (snow - melt) * 1000 / RHO_I


def simulate(bed0, climate, dx, years, dt, sea_level=-120.0):
    """Grow ice from nothing on `bed0` under LGM `climate` and return (ice thickness, depressed bed)."""
    T_ref, P, z_ref, _ = climate
    # Ice slides fast over smooth, low-lying sediment (e.g. the Laurentide's southern lobes).
    mean = gaussian_filter(bed0, 50e3 / dx)
    rough = np.sqrt(np.maximum(gaussian_filter(bed0**2, 50e3 / dx) - mean**2, 0))
    slide = SLIDE * np.where((rough < 150) & (bed0 < 800), SOFT_BED, 1.0)
    slide_x, slide_y = 0.5 * (slide[:, 1:] + slide[:, :-1]), 0.5 * (slide[1:] + slide[:-1])

    H = np.zeros_like(bed0)
    bed = bed0.copy()
    ny, nx = H.shape
    for _ in range(int(years / dt)):
        s = np.maximum(bed + H, sea_level)
        # Diffusivity on cell faces (deformation + sliding), frozen for this step.
        Dx = np.zeros((ny, nx + 1))
        Dy = np.zeros((ny + 1, nx))
        Hx, Hy = 0.5 * (H[:, 1:] + H[:, :-1]), 0.5 * (H[1:] + H[:-1])
        Dx[:, 1:-1] = GAMMA * Hx ** (N + 2) * (np.diff(s, axis=1) / dx) ** 2 + slide_x * RHO_I * G * Hx**2
        Dy[1:-1] = GAMMA * Hy ** (N + 2) * (np.diff(s, axis=0) / dx) ** 2 + slide_y * RHO_I * G * Hy**2
        r = dt / dx**2
        e, w, n, so = (r * D.ravel() for D in (Dx[:, 1:], Dx[:, :-1], Dy[1:], Dy[:-1]))

        # Implicit in H: (1 - r div D grad) H' = H + dt M + r div D grad bed
        mat = sparse.diags([1 + e + w + n + so, -e[:-1], -w[1:], -n[:-nx], -so[nx:]],
                           [0, 1, -1, nx, -nx], format="csc")
        div_bed = (Dx[:, 1:] * np.diff(bed, axis=1, append=0) - Dx[:, :-1] * np.diff(bed, axis=1, prepend=0)
                   + Dy[1:] * np.diff(bed, axis=0, append=0) - Dy[:-1] * np.diff(bed, axis=0, prepend=0))
        M = mass_balance(s, T_ref, P, z_ref)
        H = spsolve(mat, (H + dt * M + r * div_bed).ravel()).reshape(ny, nx)

        H = np.maximum(H, 0.0)
        H[bed < sea_level - MAX_WATER] = 0.0  # ice calves away over deep water
        H[[0, -1]] = H[:, [0, -1]] = 0.0
        bed += dt / TAU_ISO * (bed0 - RHO_I / RHO_M * H - bed)
    return H, bed
