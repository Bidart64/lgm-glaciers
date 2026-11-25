"""Build a standalone demo in docs/ (free map sources, no account): python -m glaciers.demo [region ...]

docs/ can be served by GitHub Pages as-is. It holds the viewer, the listed regions' computed ice (the first is
shown first; the others surround it so its box edges don't show) and the mapped outlines clipped to them.
Run the regions first (python -m glaciers <region>) and python -m glaciers.outlines.
"""

import json
import shutil
import sys

from .regions import WEB

DOCS = WEB.parent / "docs"
START = "lost_river_range"  # the region the demo opens on


def main():
    background = ["north_america"]  # coarse ice around the regions, so ice doesn't stop at their box edges
    metas = [json.loads(p.read_text()) for p in sorted((WEB / "data").glob("*.json")) if p.name != "regions.json"]
    names = sys.argv[1:] or [m["name"] for m in metas if m["name"] not in background]
    metas = [m for m in metas if m["name"] in names + background]
    data = DOCS / "data"
    shutil.rmtree(data, ignore_errors=True)  # drop files from earlier builds
    data.mkdir(parents=True)
    for page in ("index.html", "sources.html"):
        shutil.copy(WEB / page, DOCS)
    # No Mapbox token: the viewer falls back to USGS and Sentinel-2 imagery, AWS Terrain Tiles and bundled names.
    start = START if START in names else names[0]
    (DOCS / "config.js").write_text(f"window.MAPBOX_TOKEN = ''\nwindow.LGM_START = '{start}'\n"
                                    f"window.LGM_BACKGROUND = {json.dumps(background)}\n")
    for m in metas:
        shutil.copy(WEB / "data" / f"{m['name']}.png", data)
    shutil.copy(WEB / "data" / "names.geojson", data)  # peak and town names, and their fonts (python -m glaciers.names)
    shutil.copytree(WEB / "fonts", DOCS / "fonts", dirs_exist_ok=True)
    # Coarsest first, as the viewer expects (finer regions are drawn over coarser ones).
    (data / "regions.json").write_text(json.dumps(sorted(metas, key=lambda m: -m["dx"])))

    # The mapped outlines, cut down to those in or near a listed region.
    boxes = [m["bounds"] for m in metas if m["name"] in names]
    inside = lambda line: any(w - 0.2 < x < e + 0.2 and s - 0.2 < y < n + 0.2 for x, y in line[::10] for w, s, e, n in boxes)
    outlines = json.loads((WEB / "data" / "outlines.geojson").read_text())
    for f in outlines["features"]:
        f["geometry"]["coordinates"] = [l for l in f["geometry"]["coordinates"] if inside(l)]
    (data / "outlines.geojson").write_text(json.dumps(outlines, separators=(",", ":")))
    size = sum(p.stat().st_size for p in DOCS.rglob("*") if p.is_file()) / 1e6
    print(f"wrote {DOCS} ({size:.1f} MB); try it: python -m http.server -d docs")


if __name__ == "__main__":
    main()
