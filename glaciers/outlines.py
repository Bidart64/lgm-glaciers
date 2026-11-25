"""Mapped (field-evidence) LGM ice extent, to compare with the model.

Ice sheets: Batchelor et al. 2019, Nature Communications (best estimate). https://osf.io/7jen3/
Mountain glaciers: GLACIMONTIS, Lima et al. 2026, Scientific Data. https://doi.org/10.5281/zenodo.15600659

python -m glaciers.outlines  ->  web/data/outlines.geojson for the viewer
"""

import json
import pickle
import struct
import urllib.request
import zipfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

from .dem import CACHE

BATCHELOR = {"shp": "https://osf.io/download/9yhdv/", "dbf": "https://osf.io/download/gzkwc/"}
GLACIMONTIS = "https://zenodo.org/api/records/15600659/files/GLACIMONTIS_SHP.zip/content"


def _rings(shp, keep=None):
    """Polygon rings (n, 2) from a shapefile; `keep[i]` selects records."""
    pos, i, out = 100, 0, []
    while pos < len(shp):
        size = struct.unpack(">i", shp[pos + 4:pos + 8])[0] * 2
        rec = shp[pos + 8:pos + 8 + size]
        if struct.unpack("<i", rec[:4])[0] == 5 and (keep is None or keep[i]):
            nparts, npts = struct.unpack("<2i", rec[36:44])
            parts = list(struct.unpack(f"<{nparts}i", rec[44:44 + 4 * nparts])) + [npts]
            xy = np.frombuffer(rec, "<f8", 2 * npts, 44 + 4 * nparts).reshape(-1, 2)
            out += [xy[a:z] for a, z in zip(parts, parts[1:])]
        pos += 8 + size
        i += 1
    return out


def _dbf_column(dbf, name):
    n, header, length = struct.unpack("<IHH", dbf[4:12])
    off = 1
    for i in range(32, header - 1, 32):
        width = dbf[i + 16]
        if dbf[i:i + 11].split(b"\0")[0].decode() == name:
            return [dbf[header + r * length + off:header + r * length + off + width].decode("latin1").strip()
                    for r in range(n)]
        off += width


def load():
    """{source: [ring (n, 2) lon/lat]}; outer rings are clockwise, holes anticlockwise (shapefile rule). Cached."""
    path = CACHE / "outlines.pkl"
    if not path.exists():
        get = lambda url: urllib.request.urlopen(url, timeout=600).read()
        dbf = get(BATCHELOR["dbf"])
        sheets = _rings(get(BATCHELOR["shp"]), [f == "Ice" for f in _dbf_column(dbf, "FEATURE")])
        z = CACHE / "glacimontis.zip"
        if not z.exists():
            z.write_bytes(get(GLACIMONTIS))
        with zipfile.ZipFile(z) as f:
            merc = _rings(f.read(next(n for n in f.namelist() if n.endswith("FilteredReconstructedPaleoglaciers.shp"))))
        mountains = [np.degrees(np.column_stack([x / 6378137, 2 * np.arctan(np.exp(y / 6378137)) - np.pi / 2]))
                     for x, y in (r.T for r in merc)]  # Web Mercator -> degrees
        path.write_bytes(pickle.dumps({"Batchelor et al. 2019": sheets, "GLACIMONTIS 2026": mountains}))
    return pickle.loads(path.read_bytes())


def rasterize(proj, x, y, sheets=True):
    """Mapped ice on a region grid (bool, same orientation as the model arrays); `sheets=False` leaves out the
    ice-sheet outline, which is drawn too coarsely to use in some mountains."""
    img = Image.new("1", (len(x), len(y)))
    draw = ImageDraw.Draw(img)
    dx = x[1] - x[0]
    for ring in (r for source, rs in load().items() if sheets or "Batchelor" not in source for r in rs):
        px, py = proj.forward(*ring.T)
        c, r = (px - x[0]) / dx, (py - y[0]) / dx
        if len(ring) < 3 or c.max() < 0 or r.max() < 0 or c.min() > len(x) or r.min() > len(y):
            continue
        clockwise = np.sum((ring[1:, 0] - ring[:-1, 0]) * (ring[1:, 1] + ring[:-1, 1])) > 0
        draw.polygon(list(zip(c, r)), fill=int(clockwise))
    return np.array(img)


def _simplify(ring, tol):
    """Douglas-Peucker, in degrees."""
    keep = np.zeros(len(ring), bool)
    keep[[0, -1]] = True
    stack = [(0, len(ring) - 1)]
    while stack:
        a, z = stack.pop()
        if z <= a + 1:
            continue
        p, d = ring[a + 1:z] - ring[a], ring[z] - ring[a]
        n = np.hypot(*d)
        dist = np.abs(d[0] * p[:, 1] - d[1] * p[:, 0]) / n if n else np.hypot(*p.T)  # closed ring: a == z point
        k = int(dist.argmax())
        if dist[k] > tol:
            keep[a + 1 + k] = True
            stack += [(a, a + 1 + k), (a + 1 + k, z)]
    return ring[keep]


def main():
    tolerance = {"Batchelor et al. 2019": 0.01, "GLACIMONTIS 2026": 0.004}  # degrees
    features = []
    for source, rings in load().items():
        lines = [np.round(_simplify(r, tolerance[source]), 3).tolist() for r in rings if len(r) > 3]
        features.append({"type": "Feature", "properties": {"source": source},
                         "geometry": {"type": "MultiLineString", "coordinates": [l for l in lines if len(l) > 3]}})
    out = Path(__file__).resolve().parent.parent / "web" / "data" / "outlines.geojson"
    out.write_text(json.dumps({"type": "FeatureCollection", "features": features}, separators=(",", ":")))
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
