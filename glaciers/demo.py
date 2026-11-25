"""Build a standalone demo in docs/ (free map sources, no account): python -m glaciers.demo lost_river idaho

docs/ can be served by GitHub Pages as-is. It holds the viewer, the listed regions' computed ice (the first is
shown first; the others surround it so its box edges don't show) and the mapped outlines clipped to them.
Run the regions first (python -m glaciers <region>) and python -m glaciers.outlines.
"""

import json
import shutil
import sys

from .regions import WEB

DOCS = WEB.parent / "docs"


def main():
    names = sys.argv[1:] or ["lost_river", "idaho", "wind_river", "sawtooth", "rainier"]
    background = ["north_america"]  # coarse ice around the regions, so ice doesn't stop at their box edges
    metas = [json.loads((WEB / "data" / f"{r}.json").read_text()) for r in names + background]
    data = DOCS / "data"
    shutil.rmtree(data, ignore_errors=True)  # drop files from earlier builds
    data.mkdir(parents=True)
    for page in ("index.html", "sources.html"):
        shutil.copy(WEB / page, DOCS)
    # No Mapbox token: the viewer falls back to USGS imagery, AWS Terrain Tiles and OpenFreeMap names.
    (DOCS / "config.js").write_text(f"window.MAPBOX_TOKEN = ''\nwindow.LGM_START = '{names[0]}'\n"
                                    f"window.LGM_REGIONS = {json.dumps(names)}\nwindow.LGM_BACKGROUND = {json.dumps(background)}\n")
    for r in names + background:
        shutil.copy(WEB / "data" / f"{r}.png", data)
    # Coarsest first, as the viewer expects (finer regions are drawn over coarser ones).
    (data / "regions.json").write_text(json.dumps(sorted(metas, key=lambda m: -m["dx"])))

    w, s, e, n = (f(m["bounds"][i] for m in metas[:len(names)]) for i, f in enumerate((min, min, max, max)))
    inside = lambda line: any(w - 0.2 < x < e + 0.2 and s - 0.2 < y < n + 0.2 for x, y in line)
    outlines = json.loads((WEB / "data" / "outlines.geojson").read_text())
    for f in outlines["features"]:
        f["geometry"]["coordinates"] = [l for l in f["geometry"]["coordinates"] if inside(l)]
    (data / "outlines.geojson").write_text(json.dumps(outlines, separators=(",", ":")))
    size = sum(p.stat().st_size for p in DOCS.rglob("*") if p.is_file()) / 1e6
    print(f"wrote {DOCS} ({size:.1f} MB); try it: python -m http.server -d docs")


if __name__ == "__main__":
    main()
