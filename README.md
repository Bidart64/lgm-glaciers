# LGM Glaciers

Computes where glaciers should have been at the Last Glacial Maximum (~21,000 years ago) in Idaho, the
Wind River Range and around Mount Rainier, and shows them in 3D on satellite imagery.

**Live demo:** https://bidart64.github.io/lgm-glaciers/ (Idaho, the Lost River, Sawtooth and Wind River ranges, Mount Rainier;
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

The single ranges run on 200 m grids, each with its own fit (`lost_river`, `sawtooth`, `wind_river`, `rainier`):

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
6. **Export** (`glaciers/export.py`): one small PNG per region with the ice surface (smoothed over 1.5
   cells: ice cannot hold bumps much smaller than its thickness), thickness, mass balance, crevassing and
   flow lines (`web/data/<region>.png` plus `.json`).
7. **Viewer** (`web/index.html`, MapLibre): for each map tile, web workers take present-day terrain (~10 m
   in the US), raise it to the computed ice surface wherever the ice buries it, and draw that as 3D terrain
   with a painted ice overlay. The model runs at 200-500 m, but ice margins, nunataks and valley glaciers
   appear at terrain resolution: thick ice is flat across a valley, thin ice follows the fine terrain and
   stays off slopes too steep to hold it. The ice is lit from its own surface and painted from the model:
   snow above the snowline, bare blue ice below it, bands along the flow lines, crevasses where the ice
   stretches, bergschrunds under steep headwalls, debris towards the snouts, and perennial snow on ice-free ground above the snowline.
8. **Validation** (`glaciers/outlines.py`): mapped LGM extent from field evidence, from Batchelor et al. 2019
   for the ice sheets and GLACIMONTIS (Lima et al. 2026) for mountain glaciers. Each run reports the
   overlap between computed and mapped ice (intersection ÷ union), and the viewer can draw the mapped
   outlines over the model.

## Results and limitations

| Region | Grid | Computed ice | Mapped ice | Overlap | Held-out range |
|---|---|---|---|---|---|
| North America (no correction) | 20 km | 16.7M km² | 18.6M km² | 84% | – |
| Idaho | 500 m | 160,000 km² | 147,000 km² | 62% | Lemhi Range 28%, Pioneer Mountains 64% |
| Sawtooth Range | 200 m | 3,800 km² | 5,300 km² | 60% | White Clouds 68% |
| Lost River Range | 200 m | 3,400 km² | 4,000 km² | 47% | Pioneer Mountains 68% |
| Wind River Range | 200 m | 10,100 km² | 9,800 km² | 71% | East-side valleys 46% |
| Mount Rainier | 200 m | 4,800 km² | 1,700 km² (valley tongues only) | 35% below 1,100 m | Goat Rocks 68% |

For comparison, a published model of the Alps at the LGM reaches about 72% at best, and automated
reconstructions average about 48%.

- Idaho's overlap is 74% north of 47°N, where the ice-sheet lobes are, and 42% in the mountain ranges to
  the south; the 200 m runs of single ranges do better there.
- The overlap is partly a fit: the correction was chosen using the mapped outlines. The held-out ranges
  are the independent test.
- Ice thickness is not constrained by the mapped outlines. The model buries most of the Sawtooths and the
  Wind Rivers under ice more than 1 km thick. Field studies do describe a plateau ice cap with piedmont
  lobes on the Wind Rivers, and Sawtooth lobes reaching up to 20 km into Stanley Basin, but give no
  thickness to check against.
- **Mount Rainier** stands at its pre-collapse height: a cone is put back on today's summit up to about
  4,900 m (3.4 km³ of rock), in the model and in the 3D view. Around Rainier the outlines draw only the
  valley tongues: nothing on the mountain itself, which still carries glaciers, nor on slopes above the
  ice-age snowline. So the fit and the score use ground below 1,100 m only, where the model covers 73% of
  the mapped tongues but also puts ice on the ridges between them. The ice-sheet outline is not used there
  (it covers the Cascades whole), and the Puget lobe of the ice sheet, which reached the north-west corner
  a few thousand years later, is not modelled.
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

## Sources

Everything used to recreate the LGM ice, and what each is used for.

**Field evidence of LGM ice extent** (fitting and scoring; the viewer's "mapped extent" overlay)
- GLACIMONTIS, reconstructed mountain paleoglaciers worldwide: Lima et al. 2026, *Scientific Data*.
  Data: https://doi.org/10.5281/zenodo.15600659
- Northern Hemisphere ice sheets, LGM best estimate: Batchelor et al. 2019, "The configuration of Northern
  Hemisphere ice sheets through the Quaternary", *Nature Communications* 10, 3713. Data: https://osf.io/7jen3/

**Ice-age climate**
- PMIP4 / CMIP6 `lgm` and `piControl` simulations from AWI-ESM-1-1-LR, INM-CM4-8, MIROC-ES2L and
  MPI-ESM1-2-LR, downloaded from the Earth System Grid Federation (DKRZ node). Overview: Kageyama et al. 2021,
  "The PMIP4 Last Glacial Maximum experiments", *Climate of the Past* 17, 1065-1089.
- Tested and not used: CHELSA-TraCE21k (Karger et al. 2023, *Climate of the Past*).

**Present-day climate**
- PRISM 1991-2020 monthly normals, 800 m: PRISM Climate Group, Oregon State University, https://prism.oregonstate.edu
- WorldClim 2.1 (ice-sheet grid): Fick & Hijmans 2017, *International Journal of Climatology* 37, 4302-4315.

**Terrain**
- Terrain Tiles (Mapzen, hosted as AWS Open Data; 10-30 m sources in the US): https://registry.opendata.aws/terrain-tiles/
- ETOPO 2022 bedrock (ice-sheet grid): NOAA National Centers for Environmental Information,
  https://doi.org/10.25921/fd45-gt74

**Mount Rainier's lost summit**
- The Osceola collapse about 5,600 years ago: Vallance & Scott 1997, "The Osceola Mudflow from Mount Rainier",
  *GSA Bulletin* 109, 143-163. The pre-collapse height of about 4,900 m (16,000 ft) is the commonly quoted
  USGS estimate; the cone put back here is a plain cone of about that height, not a geological reconstruction.

**Checks against field studies** (not model inputs)
- Wind River ice cap: Birkel et al. 2012, "Climate inferences from a glaciological reconstruction of the late
  Pleistocene Wind River Ice Cap", *Arctic, Antarctic, and Alpine Research* 44, 265-276.
- Sawtooth Valley glaciation: Thackray et al. 2004, *Geology* 32, 225-228; USGS Open-File Report 2004-1222.

**Methods**
- Shallow-ice approximation: Hutter 1983, *Theoretical Glaciology*.
- Positive-degree-day melt with the closed-form daily-variability integral: Calov & Greve 2005,
  *Journal of Glaciology* 51, 173-175.

**Viewer**
- MapLibre GL JS; imagery from the USGS National Map; place names from OpenFreeMap / OpenMapTiles /
  OpenStreetMap; with a token, imagery, terrain and names from Mapbox instead.
