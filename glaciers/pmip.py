"""LGM climate anomalies from the PMIP4 (CMIP6 `lgm`) ensemble, applied to present climate (delta method).

Temperature anomalies are reduced to sea level with each model's own orography, so the height of the
models' LGM ice sheets is not counted twice: the ice-flow model applies its own lapse rate to its own ice.
"""

import io
import json
import urllib.parse
import urllib.request

import h5py
import numpy as np
from scipy.ndimage import map_coordinates

from .dem import CACHE

SEARCH = "https://esgf-data.dkrz.de/esg-search/search"
MODELS = {"AWI-ESM-1-1-LR": "r1i1p1f1", "INM-CM4-8": "r1i1p1f1", "MIROC-ES2L": "r1i1p1f2", "MPI-ESM1-2-LR": "r1i1p1f1"}
YEARS = 30      # climatology length (last years of each run)
LAPSE = 0.0065  # K per m, as in the ice model
GRID = np.arange(-89.5, 90), np.arange(-179.5, 180)  # common 1° grid


def _files(model, experiment, var):
    q = {"type": "File", "project": "CMIP6", "latest": "true", "distrib": "false", "limit": 1000,
         "format": "application/solr+json", "fields": "url", "source_id": model, "experiment_id": experiment,
         "member_id": MODELS[model], "variable_id": var, "table_id": "fx" if var == "orog" else "Amon"}
    docs = json.load(urllib.request.urlopen(f"{SEARCH}?{urllib.parse.urlencode(q)}", timeout=120))["response"]["docs"]
    urls = {u.split("|")[0] for d in docs for u in d["url"] if u.endswith("HTTPServer")}
    return sorted(urls, key=lambda u: u.rsplit("_", 1)[-1])  # chronological by the date range in the name


def _read(url, var):
    print("  ", url.rsplit("/", 1)[-1])
    with urllib.request.urlopen(url, timeout=600) as r, h5py.File(io.BytesIO(r.read())) as f:
        return f[var][...].astype(np.float64), f["lat"][...], f["lon"][...]


def _regrid(a, lat, lon):
    """Bilinear onto the common 1° grid (periodic in longitude)."""
    if lat[0] > lat[-1]:
        lat, a = lat[::-1], a[..., ::-1, :]
    lon = np.concatenate([lon, lon[:1] + 360])
    a = np.concatenate([a, a[..., :1]], -1)
    la, lo = np.meshgrid(*GRID, indexing="ij")
    rc = [np.interp(la, lat, np.arange(len(lat))), np.interp(lo % 360, lon, np.arange(len(lon)))]
    return np.stack([map_coordinates(m, rc, order=1) for m in a.reshape(-1, *a.shape[-2:])]).reshape(*a.shape[:-2], *la.shape)


def _climatology(model, experiment, var):
    """Mean annual cycle (12, lat, lon) over the last YEARS of a run."""
    chunks, n = [], 0
    for url in reversed(_files(model, experiment, var)):
        a, lat, lon = _read(url, var)
        chunks.insert(0, a)
        n += len(a)
        if n >= 12 * YEARS:
            break
    a = np.concatenate(chunks)[-12 * YEARS:]
    return _regrid(a.reshape(-1, 12, *a.shape[1:]).mean(0), lat, lon)


def _orog(model, experiment):
    a, lat, lon = _read(_files(model, experiment, "orog")[0], "orog")
    return _regrid(np.where(np.abs(a) > 1e10, 0, a), lat, lon)


def anomalies():
    """Ensemble-mean (sea-level temperature change K, precipitation ratio), each (12, 180, 360). Cached."""
    path = CACHE / "pmip" / "anomalies.npz"
    if not path.exists():
        dT, ratio = [], []
        for m in MODELS:
            print(f"PMIP4 {m}")
            t = {e: _climatology(m, e, "tas") + LAPSE * _orog(m, e) for e in ("lgm", "piControl")}
            p = {e: _climatology(m, e, "pr") * 86400 for e in ("lgm", "piControl")}  # mm/day
            dT.append(t["lgm"] - t["piControl"])
            ratio.append((p["lgm"] + 0.1) / (p["piControl"] + 0.1))
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, dT=np.mean(dT, 0), ratio=np.clip(np.mean(ratio, 0), 0.1, 3))
    d = np.load(path)
    return d["dT"], d["ratio"]


def sample(lon, lat):
    """Anomalies at arbitrary points: (dT, precipitation ratio), each (12, *lon.shape)."""
    rc = [lat + 89.5, (lon + 179.5) % 360]
    return tuple(np.stack([map_coordinates(m, rc, order=1, mode="grid-wrap") for m in a]) for a in anomalies())
