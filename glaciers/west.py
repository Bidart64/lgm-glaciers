"""The western United States, patch by patch: python -m glaciers.west [--jobs 7] [--memory 11] [--redo]

Only about 4% of the West carried mountain glaciers at the LGM, so instead of one grid over all of it, each group
of mapped glaciers gets its own 200 m region: a box around the group with BUFFER to spare, so that the ice can
grow beyond its outlines. The first run draws up the patches (glaciers/patches.json, named after the mountains
they hold); every run then fits and computes the patches not yet done, several at a time.

Ice-sheet lobes are not part of this: they come from the north_america run.
"""

import argparse
import io
import json
import re
import subprocess
import sys
import threading
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from PIL import Image, ImageDraw
from scipy.ndimage import binary_dilation, find_objects, label

from . import dem, names, outlines
from .regions import PATCHES, WEB

WEST, EAST, SOUTH, NORTH = -125.0, -102.0, 31.0, 49.5  # the area searched for mapped glaciers
CELL = 0.02      # degrees: the raster the glaciers are grouped on (about 2 km)
BUFFER = 10.0    # km: glaciers closer than twice this share a patch, and each box reaches this far past its ice
LOGS = dem.CACHE / "logs"
# Patches are named after the mountain range nearest the middle of their ice, which for the big ones is some
# minor local range. These get the name people know them by: a patch containing the point takes the name.
NAMED = [(-114.9, 46.5, "bitterroot_range"), (-110.5, 44.5, "yellowstone_and_wind_river"), (-119.3, 37.8, "sierra_nevada"),
         (-106.4, 39.5, "central_colorado"), (-122.0, 44.1, "oregon_cascades"), (-114.95, 44.1, "sawtooth_range"),
         (-121.76, 46.85, "mount_rainier"), (-111.7, 40.6, "wasatch_range"), (-121.1, 47.6, "central_cascades"),
         (-110.4, 40.75, "uinta_mountains"), (-123.0, 41.0, "klamath_mountains"), (-113.78, 44.14, "lost_river_range"),
         (-121.5, 40.49, "lassen_peak"), (-121.3, 48.7, "north_cascades"), (-113.8, 48.7, "glacier_national_park"),
         (-122.19, 41.41, "mount_shasta"), (-110.8, 43.74, "teton_range")]


def ranges():
    """[(lon, lat, name)] of the mountain ranges in GeoNames' United States table, to name the patches after."""
    path = dem.CACHE / "geonames" / "US.zip"
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(names.get(names.GEONAMES.format("US")))
    with zipfile.ZipFile(path) as z, io.TextIOWrapper(z.open("US.txt"), "utf-8") as rows:
        fields = (row.split("\t") for row in rows)
        return [(float(f[5]), float(f[4]), f[1]) for f in fields if f[7] == "MTS"]


def draw_patches():
    """{name: (centre lon, centre lat, half-width km, half-height km)}, largest first."""
    nx, ny = round((EAST - WEST) / CELL), round((NORTH - SOUTH) / CELL)
    img = Image.new("1", (nx, ny))
    draw = ImageDraw.Draw(img)
    for ring in outlines.load()["GLACIMONTIS 2026"]:  # mountain glaciers, rasterized as in outlines.rasterize
        lon, lat = ring.T
        if len(ring) < 3 or lon.max() < WEST or lon.min() > EAST or lat.max() < SOUTH or lat.min() > NORTH:
            continue
        clockwise = np.sum((lon[1:] - lon[:-1]) * (lat[1:] + lat[:-1])) > 0
        draw.polygon(list(zip((lon - WEST) / CELL, (NORTH - lat) / CELL)), fill=int(clockwise))
    ice = np.array(img)
    reach = round(BUFFER / (CELL * 111.2))  # the buffer, in raster cells
    groups, _ = label(binary_dilation(ice, iterations=reach))  # glaciers whose buffers touch form one group
    known, patches = ranges(), {}
    boxes = sorted(find_objects(groups), key=lambda box: -(box[0].stop - box[0].start) * (box[1].stop - box[1].start))
    for rows, cols in boxes:
        lat0 = NORTH - (rows.start + rows.stop) / 2 * CELL
        lon0 = WEST + (cols.start + cols.stop) / 2 * CELL
        hh = (rows.stop - rows.start) * CELL * 111.2 / 2
        hw = (cols.stop - cols.start) * CELL * 111.2 * np.cos(np.radians(lat0)) / 2
        # Named after the mountain range nearest the middle of its ice.
        r, c = np.nonzero(ice[rows, cols])
        lat_i, lon_i = NORTH - (rows.start + r.mean()) * CELL, WEST + (cols.start + c.mean()) * CELL
        _, title = min(((lon - lon_i) ** 2 * np.cos(np.radians(lat_i)) ** 2 + (lat - lat_i) ** 2, name) for lon, lat, name in known)
        inside = lambda lon, lat: abs(lon - lon0) * 111.2 * np.cos(np.radians(lat0)) < hw and abs(lat - lat0) * 111.2 < hh
        title = next((name for lon, lat, name in NAMED if inside(lon, lat)), title)
        name = base = re.sub(r"[^a-z]+", "_", title.lower()).strip("_")
        for n in range(2, 99):  # two patches can be nearest the same range
            if name not in patches:
                break
            name = f"{base}_{n}"
        patches[name] = (round(lon0, 3), round(lat0, 3), int(np.ceil(hw)), int(np.ceil(hh)))
    return patches


def compute(name):
    """Fit one patch's climate correction and compute its ice, logging to data/logs/<name>.log."""
    with open(LOGS / f"{name}.log", "w") as log:
        for module in ("glaciers.fit", "glaciers"):
            if subprocess.run([sys.executable, "-u", "-m", module, name], stdout=log, stderr=subprocess.STDOUT).returncode:
                return f"{name}: failed, see {log.name}"
    meta = json.loads((WEB / "data" / f"{name}.json").read_text())
    return f"{name}: ice {meta['ice_km2']:,} km2, mapped {meta['mapped_km2']:,} km2, overlap {meta['iou']:.0%}"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--jobs", type=int, default=7, help="most patches computed at once")
    p.add_argument("--memory", type=float, default=11, help="GB the patches running at once may need between them")
    p.add_argument("--redo", action="store_true", help="recompute patches that already have results")
    a = p.parse_args()
    if not PATCHES.exists():
        PATCHES.write_text(json.dumps(draw_patches(), indent=1))
        print(f"wrote {PATCHES}")
    patches = json.loads(PATCHES.read_text())
    todo = [name for name in patches if a.redo or not (WEB / "data" / f"{name}.json").exists()]
    print(f"{len(patches)} patches, {len(todo)} to compute")
    if todo:
        LOGS.mkdir(parents=True, exist_ok=True)
        dem._prism(slice(0, 1), slice(0, 1))  # fetch the climate normals once, before the jobs all ask for them
        # A patch needs about 1 GB per million cells while it runs. Each starts once the ones already running
        # leave it room (or nothing is running), so the few big patches go a couple at a time and the many small
        # ones as fast as there are cores for them.
        need = {name: (2 * hw * 5 + 1) * (2 * hh * 5 + 1) / 1e6 for name, (_, _, hw, hh) in patches.items()}  # GB
        room, used = threading.Condition(), [0.0]

        def job(name):
            with room:
                room.wait_for(lambda: used[0] == 0 or used[0] + need[name] <= a.memory)
                used[0] += need[name]
            try:
                return compute(name)
            finally:
                with room:
                    used[0] -= need[name]
                    room.notify_all()

        with ThreadPoolExecutor(a.jobs) as pool:  # threads that each wait on their own processes
            for n, done in enumerate(as_completed([pool.submit(job, name) for name in todo]), 1):
                print(f"[{n}/{len(todo)}] {done.result()}", flush=True)


if __name__ == "__main__":
    main()
