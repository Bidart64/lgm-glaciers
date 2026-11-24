# LGM Glaciers

Computes where glaciers and ice sheets should have been at the Last Glacial Maximum (~21,000 years ago)
and shows them in 3D on satellite imagery.

**Live demo:** https://bidart64.github.io/lgm-glaciers/ (Idaho, the Lost River, Sawtooth and Wind River ranges;
no account needed)

## Run

```bash
uv venv && uv pip install -e .                  # numpy, scipy, pillow, h5py
python -m glaciers alps                         # ~40 s after the first download; also: pyrenees, new_zealand,
                                                # patagonia, europe, north_america
python -m glaciers.outlines                     # mapped LGM extent: scores every run, and the viewer's overlay
cp web/config.example.js web/config.js          # optional: a Mapbox token (Mapbox imagery, terrain, names)
python -m http.server -d web                    # open http://localhost:8000
```

High-resolution Idaho and the Lost River Range (run in this order):

```bash
python -m glaciers north_america                # ~4 min: supplies the Cordilleran ice entering Idaho
python -m glaciers.fit idaho                    # ~25 min: fits Idaho's climate correction to the mapped ice
python -m glaciers idaho                        # ~17 min: 500 m grid
python -m glaciers lost_river                   # ~3 min: 200 m grid, uses Idaho's correction
python -m glaciers.fit sawtooth                 # ~35 min: Sawtooth Range correction at 200 m (Idaho's baseline)
python -m glaciers sawtooth                     # ~5 min: 200 m grid
python -m glaciers.demo                         # standalone demo in docs/ (free map sources, GitHub Pages)
```

Wind River Range, Wyoming (no parent run needed):

```bash
python -m glaciers.fit wind_river               # fits its own climate correction at 500 m
python -m glaciers wind_river                   # ~5 min: 200 m grid
```

Without a Mapbox token the viewer uses free sources: USGS National Map imagery (US only), AWS Terrain
Tiles and OpenFreeMap names. The demo in `docs/` runs that way, so anyone can open it. Controls follow
`../maps`: reset view, 3D terrain, viewing angle, place names and mapped extent. Tilt and rotate with
Ctrl+drag, Option+drag (Alt on Windows) or right-drag. The URL keeps the view, so views can be bookmarked and shared.

Options: `--dT -1` (extra cooling on top of the LGM climate), `--sea-level -130`, `--years 5000`. The
first runs download about 2.5 GB (PMIP4 climate-model output, PRISM normals, glacier outlines, bedrock)
and cache about 400 MB in `data/`.

## How it works

1. **Bedrock**: ETOPO 2022 bedrock, 1 arc-minute (30 arc-seconds for grids finer than 2 km), with
   today's Greenland and Antarctic ice removed. It is fetched from NOAA over OPeNDAP and resampled onto
   an equal-area grid per region. A nested region (Idaho, inside North America) holds its edges at the
   parent run's ice thickness, so ice sheets can flow in from outside the box.
2. **Climate**: present-day monthly temperature and precipitation from WorldClim 2.1 (≈18 km), plus the
   LGM-minus-preindustrial change from the four PMIP4 models with monthly LGM output on ESGF
   (AWI-ESM-1-1-LR, INM-CM4-8, MIROC-ES2L, MPI-ESM1-2-LR; 30-year climatologies, ensemble mean). This is
   the standard "delta method". Temperature changes are reduced to sea level with each model's own LGM
   orography, so the climate models' ice-sheet height isn't counted twice, and precipitation changes are
   applied as ratios. Temperatures are then corrected to the computed ice surface with a lapse rate
   (`glaciers/pmip.py`).
3. **Mass balance**: a positive-degree-day model computes snowfall minus melt, with 60% of snowmelt
   refreezing.
4. **Ice flow**: a shallow-ice model (Glen's flow law plus basal sliding, five times faster on smooth
   lowland sediment), solved semi-implicitly. Ice grows from nothing until it settles. Ice can ground in
   seas shallower than 500 m, and the bed sinks under the ice load.
5. **Export**: one small PNG per region holding the ice surface elevation and the exposed shelf
   (`web/data/<region>.png` plus `.json`).
6. **Viewer** (`web/index.html`, MapLibre): for each map tile, the browser takes present-day terrain
   (Mapbox Terrain-DEM: ~10 m in the US, 30 m most elsewhere), raises it to the computed ice surface
   wherever the ice buries it, and draws that as 3D terrain with an ice overlay. The model runs at 1–20 km, but ice margins, nunataks and valley
   glaciers appear at terrain resolution. This is the standard downscaling trick: ice surfaces are
   smooth, while the ground underneath is not. Within 3 km of a margin, ice thins as `√(d / 3 km)`, where
   d is the distance to the margin, like an ice-sheet edge profile, so margins taper instead of ending in
   cliffs. Ice ending against higher rock keeps its full thickness.
7. **Validation** (`glaciers/outlines.py`): mapped LGM extent from field evidence, from Batchelor et al.
   2019 for the Northern Hemisphere ice sheets and GLACIMONTIS (Lima et al. 2026) for mountain glaciers
   worldwide. Each run reports the overlap between computed and mapped ice (intersection ÷ union), and
   the viewer can draw the mapped outlines over the model.

8. **Climate correction (Idaho)** (`glaciers/fit.py`): the climate models' LGM summer cooling has the
   wrong regional shape over Idaho, about −6.5 °C in the central ranges (glacier modelling in the northern
   Rockies implies about −8.5 to −10 °C) and −9 to −11 °C near the ice sheets. CHELSA-TraCE21k, a 1 km
   downscaled alternative, is no better (−5.4 to −6 °C centrally). So a smooth temperature correction on a
   1° grid is fitted to the mapped ice, like climate inversions from moraines. The Lemhi and Pioneer
   ranges are held out of the fit to test it. Idaho's ice is therefore a reconstruction constrained by
   field evidence, not a pure prediction.

Code: `glaciers/dem.py` (data), `glaciers/pmip.py` (LGM climate), `glaciers/model.py` (physics),
`glaciers/outlines.py` (validation), `glaciers/fit.py` (correction), `glaciers/regions.py`, `glaciers/export.py`,
`glaciers/__main__.py` (CLI), `web/index.html`.

## Calibration and limitations

| Region | Extra cooling | Computed ice | Mapped ice | Overlap |
|---|---|---|---|---|
| North America | 0 | 16.3M km² | 18.6M km² | 83% |
| Europe | 0 | 6.4M km² | 9.2M km² | 67% |
| Alps (incl. Jura, Black Forest, Massif Central) | −3.5 °C | 152,000 km² | 124,000 km² | 54% |
| Pyrenees | −3.5 °C | 8,900 km² | 11,600 km² | 40% |
| New Zealand | −3.5 °C | 23,000 km² | 63,000 km² | 31% |
| Idaho (500 m, PRISM climate, fitted correction) | fitted | 164,000 km² | 147,000 km² | 59% |
| Lost River Range (200 m) | fitted (Idaho's) | 4,100 km² | 4,000 km² | 38% |
| Patagonia | −3.5 °C | 589,000 km² | not in either dataset (≈ 480,000 in PATICE) | – |

The ice sheets use the PMIP4 climate unchanged. Their southern margin is at 39–42°N east of the Rockies,
interior Alaska is ice-free apart from the ranges, and Siberia is ice-free. The 2–3° climate models
can't resolve mountain ranges and underestimate their LGM cooling. One shared offset for the mountain
regions maximizes their mean overlap (34% at −2.5, 42% at −3.5, 43% at −4, where the Alps grow to 1.7×
the mapped area). The same offset everywhere makes North America too large (18.3M km² at −2.5). Because
that offset was chosen using the mapped outlines, the mountain overlaps are a fit, not an independent
test; the ice-sheet overlaps are independent.

The Patagonian reconstruction (PATICE, Davies et al. 2020) is on Mendeley Data, which needs a login to
download, so it isn't included.

Known gaps:

- The marine Barents–Kara ice sheet barely forms (7% of the Barents Sea), which is most of Europe's
  missing overlap. Its growth depends on
  ice-shelf and grounding-line processes that a shallow-ice model doesn't include.
- Within mountain regions, ice is in the wrong places as well as the wrong amount: New Zealand is too
  small and the Alps too large at the same offset, which points to regional climate biases. The offset
  also puts thin ice caps on German uplands around 50°N, which were not glaciated.
- The ice-sheet climate is partly circular: the PMIP4 models were run with prescribed LGM ice sheets.
- Present-day fields fill in the sea floor exposed by the lower sea level.
