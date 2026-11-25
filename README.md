# LGM Glaciers

Computes where glaciers should have been at the Last Glacial Maximum (~21,000 years ago) across the
western United States and Alaska, and shows them in 3D on satellite imagery.

**Live demo:** https://bidart64.github.io/lgm-glaciers/ (95 mountain areas from the Sierra Nevada to Glacier National Park, plus Alaska;
no account needed)

## Run

```bash
uv venv && uv pip install -e .                  # numpy, scipy, pillow, h5py
python -m glaciers.outlines                     # mapped LGM extent: scores every run, and the viewer's overlay
python -m glaciers north_america                # ~5 min: the ice sheets (the map's background, and Alaska's inflow)
python -m glaciers.west                         # ~6 h on 8 cores: every glaciated range of the western US, 200 m grids
python -m glaciers.names                        # peak and town names and fonts for the free viewer
cp web/config.example.js web/config.js          # optional: a Mapbox token (Mapbox imagery, terrain, names)
python -m http.server -d web                    # open http://localhost:8000
```

`glaciers.west` draws a box around each group of mapped glaciers (`glaciers/patches.json`, 95 patches named after
their mountains), then fits and computes each one, several at a time within a memory budget. Only about 4% of the
West carried mountain glaciers, so this costs far less than one grid over all of it. A patch can also be run alone,
and Alaska is one large region:

```bash
python -m glaciers.fit sawtooth_range           # ~1 h: fits the climate correction to the mapped ice
python -m glaciers sawtooth_range               # ~5 min
python -m glaciers sawtooth_range --mapped      # instead: fill the mapped extent with ice (see below)
python -m glaciers.fit alaska && python -m glaciers alaska   # 1.5 km grid, joined to the ice sheet
python -m glaciers.demo                         # standalone demo in docs/ (free map sources, GitHub Pages)
```

Without a Mapbox token the viewer uses free sources: USGS National Map imagery (US only), AWS Terrain
Tiles and GeoNames place names. The demo in `docs/` runs that way, so anyone can open it. Controls follow
`../maps`: reset view, 3D terrain, glaciers on/off, viewing angle, place names, region details and mapped
extent. Click the map for a point's coordinates, its elevation with the ice, and the ice thickness there. Tilt and
rotate with Ctrl+drag, Option+drag (Alt on Windows) or right-drag. The URL keeps the view, so views can be
bookmarked and shared.

The first runs download several GB (PMIP4 climate-model output, PRISM normals, glacier outlines, terrain)
and cache them in `data/` (PRISM alone is 2 GB).

## How it works

1. **Bedrock** (`glaciers/dem.py`): today's terrain from AWS Terrain Tiles (10-30 m sources), resampled onto
   an equal-area grid per region; ETOPO 2022 bedrock for the 20 km ice-sheet grid. Alaska holds its
   edges at the ice-sheet run's thickness, so the Cordilleran ice can flow in. The western patches have no
   inflow: the ice-sheet lobes that reached northern Washington, Idaho and Montana are shown only by the
   coarse North America run.
2. **Climate** (`glaciers/dem.py`, `glaciers/pmip.py`): present-day monthly temperature and precipitation
   (PRISM 800 m normals in the contiguous US; WorldClim at 4.6 km in Canada and Alaska and at 18 km for the ice sheets) plus the LGM-minus-preindustrial change from four
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
   Only ice thicker than 20 m counts as glacier; thinner ice is a film on slopes above the snowline. Part
   of each region is held out of the fit to test it: one in five quarter-degree squares, or a named range
   in Alaska. Ranges with under 500 km² of mapped ice are too small for that and hold nothing out. The ice is therefore a reconstruction constrained
   by field evidence, not a pure prediction.
6. **Export** (`glaciers/export.py`): one small PNG per region with the ice surface (smoothed over 1.5
   cells: ice cannot hold bumps much smaller than its thickness), thickness, mass balance, crevassing and
   flow lines (`web/data/<region>.png` plus `.json`).
7. **Viewer** (`web/index.html`, MapLibre): for each map tile, web workers take present-day terrain (~10 m
   in the US), raise it to the computed ice surface wherever the ice buries it, and draw that as 3D terrain
   with a painted ice overlay. The model runs at 200 m, but ice margins, nunataks and valley glaciers
   appear at terrain resolution: thick ice is flat across a valley, thin ice follows the fine terrain and
   stays off slopes too steep to hold it. The ice is lit from its own surface and painted from the model:
   snow above the snowline, bare blue ice below it, bands along the flow lines, crevasses where the ice
   stretches, bergschrunds under steep headwalls, debris towards the snouts, and perennial snow on ice-free ground above the snowline.
8. **Validation** (`glaciers/outlines.py`): mapped LGM extent from field evidence, from Batchelor et al. 2019
   for the ice sheets and GLACIMONTIS (Lima et al. 2026) for mountain glaciers. Each run reports the
   overlap between computed and mapped ice (intersection ÷ union), and the viewer can draw the mapped
   outlines over the model.

## Results and limitations

| Region | Grid | Computed ice | Mapped ice | Overlap | Held-out part |
|---|---|---|---|---|---|
| North America (no correction) | 20 km | 16.7M km² | 18.6M km² | 84% | – |
| Alaska | 1.5 km | 642,000 km² | 657,000 km² | 64% | Central Brooks Range 15% |
| Western US, 95 patches together | 200 m | 144,800 km² | 138,200 km² | 52% (weighted by mapped ice) | – |
| Yellowstone and Wind River | 200 m | 43,700 km² | 39,400 km² | 69% | 71% |
| Sierra Nevada | 200 m | 14,400 km² | 15,300 km² | 66% | 71% |
| Sawtooth Range | 200 m | 9,700 km² | 10,300 km² | 60% | 57% |
| Uinta Mountains | 200 m | 3,600 km² | 3,200 km² | 59% | 63% |
| San Juan Mountains | 200 m | 9,800 km² | 9,500 km² | 57% | 48% |
| Central Colorado | 200 m | 8,100 km² | 8,200 km² | 52% | 50% |
| Oregon Cascades | 200 m | 6,900 km² | 7,900 km² | 45% | 39% |
| Lost River Range | 200 m | 3,200 km² | 3,600 km² | 45% | 29% |
| Olympic Mountains | 200 m | 2,000 km² | 1,200 km² (valley tongues only) | 33% below 1,100 m | 45% |
| Mount Rainier | 200 m | 4,800 km² | 2,200 km² (valley tongues only) | 29% below 1,100 m | 22% |
| Bitterroot Range | 200 m | 17,800 km² | 14,400 km² | 23% | 21% |
| Central Cascades | 200 m | 3,300 km² | 1,900 km² (valley tongues only) | 21% below 1,100 m | 34% |

Overlap is intersection ÷ union of computed and mapped ice; the held-out part is one in five quarter-degree
squares, left out of the fit. Every patch's numbers are in `web/data/<patch>.json` and in the viewer's panel.

- The overlap is partly a fit: the correction was chosen using the mapped outlines. The held-out squares
  are the independent test.
- **The evidence differs from place to place**, and the scores with it:
  - In Washington's Cascades and Olympics the outlines draw only the valley tongues, not the high ground
    they flowed from, so fits and scores there use ground below 1,100 m only. The model also puts ice on
    the ridges between the tongues.
  - The ice-sheet outline counts as evidence only south of 46.5°N and east of 120°W, where it shows the
    mountain ice of central Idaho and Yellowstone. Further west it covers the Cascades whole.
  - Many small ranges are drawn as rough blobs marking that ice existed, not as glacier outlines. The
    model's glaciers can sit on the high ground beside a blob and score 0%. For the 68 patches with under
    500 km² of mapped ice (8,000 km² in all, against 5,900 km² computed), compare areas, not overlap.
- **No ice sheet in the patches.** North of about 46.5°N the Cordilleran ice sheet flowed in from Canada;
  the patches have no inflow, so the Bitterroot patch and those near the border show local glaciers only
  (hence its 23%), and nothing is patched north of 48.3°N, where the mountain-glacier outlines stop. The
  coarse North America run shows the ice sheet there.
- **Patch edges.** Each box reaches 10 km past its mapped ice. Where the computed ice grows further, it is
  cut off along the edge of the box in a straight line (north of Mount Rainier, for instance).
- About a dozen small, dry ranges (mostly in the Great Basin) reach the fit's limit of 4 °C of extra
  cooling and still grow little or no ice: 10 patches end with under 1 km².
- Ice thickness is not constrained by the mapped outlines. The model buries most of the Sawtooths and the
  Wind Rivers under ice more than 1 km thick. Field studies do describe a plateau ice cap with piedmont
  lobes on the Wind Rivers, and Sawtooth lobes reaching up to 20 km into Stanley Basin, but give no
  thickness to check against.
- **Mount Rainier** stands at its pre-collapse height: a cone is put back on today's summit up to about
  4,900 m (3.4 km³ of rock), in the model and in the 3D view. The Puget lobe of the ice sheet, which
  reached the area a few thousand years later, is not modelled.
- Where a basin is close to filling with ice the result is sensitive: a few tenths of a degree turn valley
  glaciers into an ice cap over the whole area. The fit then swings between the two; steps that reverse are
  halved, which did not settle the Central Cascades (repeated fits ended anywhere between 2,400 and
  8,200 km²), the least reliable of the large patches.
- Alaska's fit had not settled when it was stopped, and its held-out Brooks Range is poorly matched.
- **Mapped mode** (`--mapped`) gives up on prediction: ice is 4 °C colder inside the mapped extent and
  cannot survive outside it, so it fills the outlines (about 90% overlap in the ranges tried) and only its
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
- WorldClim 2.1 (outside the contiguous US: 4.6 km for mountain regions, 18 km for the ice-sheet grid): Fick & Hijmans 2017, *International Journal of Climatology* 37, 4302-4315.

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
- MapLibre GL JS; imagery from the USGS National Map (US) and Sentinel-2 cloudless by EOX (elsewhere;
  contains modified Copernicus Sentinel data 2016); peak and town names from GeoNames (geonames.org, CC BY 4.0),
  bundled with the site by `python -m glaciers.names`; with a token, imagery, terrain and names from Mapbox instead.
