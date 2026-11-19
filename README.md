# LGM Glaciers

Computes where glaciers and ice sheets should have been at the Last Glacial Maximum (~21,000 years ago)
and shows them in 3D on satellite imagery.

## Run

```bash
uv venv && uv pip install numpy scipy pillow h5py   # or: pip install -e .
python -m glaciers alps                         # ~40 s after the first download; also: pyrenees, new_zealand, patagonia, europe, north_america
python -m glaciers.outlines                     # optional: mapped LGM extent for the viewer's overlay
cp web/config.example.js web/config.js          # put your Mapbox token in it (satellite imagery)
python -m http.server -d web                    # open http://localhost:8000
```

Every computed region is shown at once. Options: `--dT -1` (extra cooling on top of the LGM climate),
`--sea-level -130`, `--years 5000`. The first run downloads about 650 MB (mostly PMIP4 climate-model
output, plus 130 MB of mapped glacier outlines) and caches about 250 MB in `data/`. The URL keeps the
current view, so views can be bookmarked and shared.

## How it works

1. **Bedrock**: ETOPO 2022 bedrock, 1 arc-minute, with today's Greenland and Antarctic ice removed.
   It is fetched from NOAA over OPeNDAP and resampled onto an equal-area grid per region.
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
6. **Viewer** (`web/index.html`, MapLibre): for each map tile, the browser takes 30 m present-day terrain
   (AWS Terrain Tiles), raises it to the computed ice surface wherever the ice buries it, and draws that
   as 3D terrain with an ice overlay. The model runs at 2–20 km, but ice margins, nunataks and valley
   glaciers appear at terrain resolution. This is the standard downscaling trick: ice surfaces are
   smooth, while the ground underneath is not. Near a margin, ice is limited to the thickness of a
   perfectly plastic ice sheet, `H = √(2τd/ρg)` with yield stress τ = 100 kPa and d the distance to the
   margin, so margins taper instead of ending in cliffs.
7. **Validation** (`glaciers/outlines.py`): mapped LGM extent from field evidence, from Batchelor et al.
   2019 for the Northern Hemisphere ice sheets and GLACIMONTIS (Lima et al. 2026) for mountain glaciers
   worldwide. Each run reports the overlap between computed and mapped ice (intersection ÷ union), and
   the viewer can draw the mapped outlines over the model.

Code: `glaciers/dem.py` (data), `glaciers/pmip.py` (LGM climate), `glaciers/model.py` (physics),
`glaciers/outlines.py` (validation), `glaciers/export.py`,
`glaciers/__main__.py` (regions and CLI), `web/index.html`.

## Calibration and limitations

| Region | Extra cooling | Computed ice | Mapped ice | Overlap |
|---|---|---|---|---|
| North America | 0 | 16.3M km² | 18.6M km² | 83% |
| Europe | 0 | 6.4M km² | 9.2M km² | 67% |
| Alps (incl. Jura, Black Forest, Massif Central) | −3.5 °C | 152,000 km² | 124,000 km² | 54% |
| Pyrenees | −3.5 °C | 8,900 km² | 11,600 km² | 40% |
| New Zealand | −3.5 °C | 23,000 km² | 63,000 km² | 31% |
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
