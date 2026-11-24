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
python -m glaciers.fit idaho                    # ~2 h: fits Idaho's climate correction to the mapped ice
python -m glaciers idaho                        # ~35 min: 500 m grid
cp web/config.example.js web/config.js          # optional: a Mapbox token (Mapbox imagery, terrain, names)
python -m http.server -d web                    # open http://localhost:8000
```

The single ranges run on 200 m grids, each with its own fit (`lost_river`, `sawtooth`, `wind_river`):

```bash
python -m glaciers.fit sawtooth                 # ~15 min
python -m glaciers sawtooth                     # ~2 min
python -m glaciers sawtooth --mapped            # instead: fill the mapped extent with ice (see below)
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
4. **Ice flow** (`glaciers/model.py`): a shallow-ice model (Glen's flow law plus slow basal sliding; with
   fast sliding, ice drains off slopes and leaves only a film). Each time
   step is one Newton step of the fully implicit equation, so both parts of the flow (diffusion of the
   surface, and thickness waves travelling with the ice) are implicit. The time step is the grid spacing
   ÷ 80 m/yr, below which results stop changing (2.5 years on a 200 m grid). Ice first settles on grids
   4× and 2× coarser, then adjusts on the full grid, which gives the same ice (within 3%) about 4× sooner. The bed
   sinks under the ice load.
5. **Climate correction** (`glaciers/fit.py`): the climate models' LGM cooling has the wrong regional shape
   over the Rockies: about −6.5 °C in central Idaho, where glacier modelling implies −8.5 to −10 °C, and
   too cold near the ice sheets. So a smooth temperature correction on a 0.5-1° grid is fitted to the mapped
   ice, like climate inversions from moraines. Every iteration grows the ice from nothing, as a run does.
   Only ice thicker than 20 m counts as glacier; thinner ice is a film on slopes above the snowline. One
   range per region is held out of the fit to test it. The ice is therefore a reconstruction constrained
   by field evidence, not a pure prediction.
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
| North America (no correction) | 20 km | 16.7M km² | 18.6M km² | 84% | – |
| Idaho | 500 m | 165,000 km² | 147,000 km² | 61% | Lemhi Range 28%, Pioneer Mountains 64% |
| Sawtooth Range | 200 m | 3,800 km² | 5,300 km² | 60% | White Clouds 68% |
| Lost River Range | 200 m | 3,400 km² | 4,000 km² | 47% | Pioneer Mountains 68% |
| Wind River Range | 200 m | 10,100 km² | 9,800 km² | 71% | East-side valleys 46% |

For comparison, a published model of the Alps at the LGM reaches about 72% at best, and automated
reconstructions average about 48%.

- Idaho's overlap is 73% north of 47°N, where the ice-sheet lobes are, and 42% in the mountain ranges to
  the south; the 200 m runs of single ranges do better there.
- The overlap is partly a fit: the correction was chosen using the mapped outlines. The held-out ranges
  are the independent test.
- Ice thickness is not constrained by the mapped outlines. The outlines show an icefield over most of the
  Sawtooth area, and the model can only cover Stanley Basin with a dome over 1 km thick that buries most
  of the range. Either the outlines are generous there or the model's ice is too stiff.
- Where a basin is close to filling with ice the result is sensitive: in the Sawtooths a few tenths of a
  degree, or 10,000 years instead of 3,000, turn valley glaciers into an ice cap over the whole area.
- **Mapped mode** (`--mapped`) gives up on prediction: ice is 4 °C colder inside the mapped extent and
  cannot survive outside it, so it fills the outlines (87-95% overlap for the three ranges) and only its
  thickness and flow are computed. It is thicker still, for the same reason.
- The ice-sheet climate is partly circular: the PMIP4 models were run with prescribed LGM ice sheets.
- Earlier versions of this model held the thickness-wave part of the flow fixed during each time step.
  That is unstable at practical time steps and created ice (in the Sawtooths at 200 m, nearly twice the
  area and over ten times the volume), which
  the fitted corrections then hid. The numbers above come from the converged solver.
