# LGM Glaciers

Computes where glaciers should have been at the Last Glacial Maximum (~21,000 years ago) in Idaho and the
Wind River Range, and shows them in 3D on satellite imagery.

**Live demo:** https://bidart64.github.io/lgm-glaciers/ (Idaho, the Lost River, Sawtooth and Wind River ranges;
no account needed)

## Run

```bash
uv venv && uv pip install -e .                  # numpy, scipy, pillow, h5py
python -m glaciers.outlines                     # mapped LGM extent: scores every run, and the viewer's overlay
python -m glaciers north_america                # ~5 min: the ice sheets, which supply the ice entering Idaho
python -m glaciers.fit idaho                    # ~40 min: fits Idaho's climate correction to the mapped ice
python -m glaciers idaho                        # ~20 min: 500 m grid
cp web/config.example.js web/config.js          # optional: a Mapbox token (Mapbox imagery, terrain, names)
python -m http.server -d web                    # open http://localhost:8000
```

The single ranges run on 200 m grids, each with its own fit (`lost_river`, `sawtooth`, `wind_river`):

```bash
python -m glaciers.fit sawtooth                 # ~5 min
python -m glaciers sawtooth                     # ~1.5 min
python -m glaciers.demo                         # standalone demo in docs/ (free map sources, GitHub Pages)
```

Without a Mapbox token the viewer uses free sources: USGS National Map imagery (US only), AWS Terrain
Tiles and OpenFreeMap names. The demo in `docs/` runs that way, so anyone can open it. Controls follow
`../maps`: reset view, 3D terrain, glaciers on/off, viewing angle, place names and mapped extent. Tilt and
rotate with Ctrl+drag, Option+drag (Alt on Windows) or right-drag. The URL keeps the view, so views can be
bookmarked and shared.

The first runs download about 2.5 GB (PMIP4 climate-model output, PRISM normals, glacier outlines, terrain)
and cache about 400 MB in `data/`.

## How it works

1. **Bedrock** (`glaciers/dem.py`): today's terrain from AWS Terrain Tiles (10-30 m sources), resampled onto
   an equal-area grid per region; ETOPO 2022 bedrock for the 20 km ice-sheet grid. Idaho holds its northern
   edge at the ice-sheet run's thickness, so the Cordilleran lobes can flow in.
2. **Climate** (`glaciers/dem.py`, `glaciers/pmip.py`): present-day monthly temperature and precipitation
   (PRISM 800 m normals; WorldClim for the ice sheets) plus the LGM-minus-preindustrial change from four
   PMIP4 climate models (AWI-ESM-1-1-LR, INM-CM4-8, MIROC-ES2L, MPI-ESM1-2-LR; ensemble mean). This is the
   standard "delta method". Temperature changes are reduced to sea level with each model's own LGM
   orography, and precipitation changes are applied as ratios. Temperatures are then corrected to the
   computed ice surface with a lapse rate.
3. **Mass balance** (`glaciers/model.py`): a positive-degree-day model computes snowfall minus melt, with
   60% of snowmelt refreezing.
4. **Ice flow** (`glaciers/model.py`): a shallow-ice model (Glen's flow law plus basal sliding). Each time
   step is one Newton step of the fully implicit equation, so both parts of the flow (diffusion of the
   surface, and thickness waves travelling with the ice) are implicit. The time step is the grid spacing
   ÷ 80 m/yr, below which results stop changing (2.5 years on a 200 m grid). Ice first settles on grids
   4× and 2× coarser, then adjusts on the full grid, which gives the same ice (within 3%) about 4× sooner. The bed
   sinks under the ice load.
5. **Climate correction** (`glaciers/fit.py`): the climate models' LGM cooling has the wrong regional shape
   over the Rockies: about −6.5 °C in central Idaho, where glacier modelling implies −8.5 to −10 °C, and
   too cold near the ice sheets. So a smooth temperature correction on a 0.5-1° grid is fitted to the mapped
   ice, like climate inversions from moraines. One range per region is held out of the fit to test it. The
   ice is therefore a reconstruction constrained by field evidence, not a pure prediction.
6. **Export** (`glaciers/export.py`): one small PNG per region with the ice surface, thickness, mass balance,
   crevassing and flow lines (`web/data/<region>.png` plus `.json`).
7. **Viewer** (`web/index.html`, MapLibre): for each map tile, web workers take present-day terrain (~10 m
   in the US), raise it to the computed ice surface wherever the ice buries it, and draw that as 3D terrain
   with a painted ice overlay. The model runs at 200-500 m, but ice margins, nunataks and valley glaciers
   appear at terrain resolution: thick ice is flat across a valley, thin ice follows the fine terrain and
   stays off slopes too steep to hold it. The ice is lit from its own surface and painted from the model:
   snow above the snowline, bare blue ice below it, bands along the flow lines, crevasses where the ice
   stretches, debris towards the snouts, and perennial snow on ice-free ground above the snowline.
8. **Validation** (`glaciers/outlines.py`): mapped LGM extent from field evidence, from Batchelor et al. 2019
   for the ice sheets and GLACIMONTIS (Lima et al. 2026) for mountain glaciers. Each run reports the
   overlap between computed and mapped ice (intersection ÷ union), and the viewer can draw the mapped
   outlines over the model.

## Results and limitations

| Region | Grid | Computed ice | Mapped ice | Overlap | Held-out range |
|---|---|---|---|---|---|
| North America (no correction) | 20 km | 16.4M km² | 18.6M km² | 84% | – |
| Idaho | 500 m | 138,000 km² | 147,000 km² | 56% | Lemhi Range 31%, Pioneer Mountains 61% |
| Sawtooth Range | 200 m | 4,900 km² | 5,300 km² | 64% | White Clouds 73% |
| Lost River Range | 200 m | 3,200 km² | 4,000 km² | 48% | Pioneer Mountains 63% |
| Wind River Range | 200 m | 9,600 km² | 9,800 km² | 73% | East-side valleys 48% |

For comparison, a published model of the Alps at the LGM reaches about 72% at best, and automated
reconstructions average about 48%.

- Idaho's overlap is 65% north of 47°N, where the ice-sheet lobes are, and 41% in the mountain ranges to
  the south; the 200 m runs of single ranges do better there.
- The overlap is partly a fit: the correction was chosen using the mapped outlines. The held-out ranges
  are the independent test.
- In the Lost River Range the ice sits too high: too much above about 2,600 m and too little below
  2,200 m, so the longest valley glaciers stop short of their mapped moraines. Changing the melt factors
  does not fix this, which points to the precipitation pattern.
- Ice thickness is not constrained by the mapped outlines. Where the fit glaciates a basin (the Sawtooth
  Valley), the ice there is probably too thick.
- The ice-sheet climate is partly circular: the PMIP4 models were run with prescribed LGM ice sheets.
- Earlier versions of this model held the thickness-wave part of the flow fixed during each time step.
  That is unstable at practical time steps and created ice (in the Sawtooths at 200 m, nearly twice the
  area and over ten times the volume), which
  the fitted corrections then hid. The numbers above come from the converged solver.
