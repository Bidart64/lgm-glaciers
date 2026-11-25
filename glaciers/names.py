"""Peak and town names for the free viewer, from GeoNames: python -m glaciers.names [region ...]

Writes web/data/names.geojson and the label fonts in web/fonts/. The viewer loads both from its own site rather
than from a map service, which content blockers may block. Run the regions first.
"""

import io
import json
import sys
import urllib.parse
import urllib.request
import zipfile

from .dem import CACHE
from .regions import REGIONS, WEB

GEONAMES = "https://download.geonames.org/export/dump/{}.zip"  # one tab-separated table per country (CC BY 4.0)
COUNTRIES = "US", "CA"
TOWNS = (50000, "city"), (2000, "town"), (50, "village")  # smallest population of each kind
FONTS = "https://tiles.openfreemap.org/fonts/{}/{}.pbf"  # Noto Sans (SIL Open Font License)
STACKS = "Noto Sans Regular", "Noto Sans Italic"
RANGES = "0-255", "256-511", "9472-9727"  # Latin, accented Latin, and the ▲ before peak names


def get(url):
    """Download; both servers turn away requests that don't say what they are."""
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "lgm-glaciers"}), timeout=600).read()


def main():
    boxes = [json.loads((WEB / "data" / f"{r}.json").read_text())["bounds"]
             for r in sys.argv[1:] or [r for r in REGIONS if r != "north_america"]]
    features = []
    for country in COUNTRIES:
        path = CACHE / "geonames" / f"{country}.zip"
        if not path.exists():
            path.parent.mkdir(parents=True, exist_ok=True)
            print(f"downloading GeoNames {country}")
            path.write_bytes(get(GEONAMES.format(country)))
        with zipfile.ZipFile(path) as z, io.TextIOWrapper(z.open(f"{country}.txt"), "utf-8") as rows:
            for row in rows:
                f = row.rstrip("\n").split("\t")  # 1 name, 4 lat, 5 lon, 6 class, 7 code, 14 population, 15 elevation, 16 terrain model
                lat, lon = float(f[4]), float(f[5])
                if not any(w < lon < e and s < lat < n for w, s, e, n in boxes):
                    continue
                if f[6] == "T" and f[7] in ("MT", "PK", "VLC"):
                    props = {"kind": "peak", "name": f[1], "ele": int(f[15] or f[16] or 0)}
                elif f[6] == "P" and int(f[14] or 0) >= TOWNS[-1][0]:
                    props = {"kind": next(kind for least, kind in TOWNS if int(f[14]) >= least), "name": f[1]}
                else:
                    continue
                features.append({"type": "Feature", "properties": props,
                                 "geometry": {"type": "Point", "coordinates": [round(lon, 5), round(lat, 5)]}})
    out = WEB / "data" / "names.geojson"
    out.write_text(json.dumps({"type": "FeatureCollection", "features": features}, ensure_ascii=False, separators=(",", ":")))
    print(f"wrote {out}: {len(features):,} names, {out.stat().st_size / 1e6:.1f} MB")
    for stack in STACKS:
        (WEB / "fonts" / stack).mkdir(parents=True, exist_ok=True)
        for chars in RANGES:
            (WEB / "fonts" / stack / f"{chars}.pbf").write_bytes(get(FONTS.format(urllib.parse.quote(stack), chars)))
    print(f"wrote {WEB / 'fonts'}")


if __name__ == "__main__":
    main()
