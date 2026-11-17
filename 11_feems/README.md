# 11 — FEEMS: fast effective migration surfaces

Estimate where gene flow is lower or higher than expected with FEEMS, a fast
successor of EEMS ([09](../09_eems/)) that also chooses the smoothing level
by cross-validation.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/make_feems_grid.R`](scripts/make_feems_grid.R) | Makes a triangular grid for the study area, finer than the grids that come with FEEMS |
| [`scripts/run_feems.py`](scripts/run_feems.py) | Reads the VCF and coordinates, builds the graph, runs the cross-validation, fits FEEMS and plots the migration surface |

**Requirements:** FEEMS (Python), R with `dggridR` and `sf` (for the grid). A
filtered VCF **without the outgroup** from [04_filtering](../04_filtering/),
the coordinates of each individual, and a habitat polygon.

---

## 1. Background

Like EEMS, FEEMS places a grid of nodes (demes) over the study area and
estimates an **effective migration** weight for each edge between nodes.
Low weights (orange) mark barriers to gene flow; high weights (blue) mark
corridors.

Differences from EEMS:

- **Fast**: it fits the model by optimisation, not MCMC (seconds to minutes
  instead of hours).
- **Smoothing (λ)**: a penalty controls how smooth the surface is. Small λ
  gives a patchy surface (risk of over-fitting); large λ a flat one. FEEMS
  chooses λ by **cross-validation**.
- **The grid** comes from a shapefile of triangles. The grids distributed
  with FEEMS are built for continental scales.

## 2. Install FEEMS

FEEMS has many dependencies. Install it in its own conda environment:

```bash
conda create -n feems -c conda-forge -c bioconda python=3.10 feems
conda activate feems
python -c "import feems; print('ok')"
```

If conda cannot solve the environment, install the dependencies first and
then FEEMS:

```bash
conda create -n feems -c conda-forge python=3.10 numpy scipy scikit-learn \
    matplotlib pyproj networkx shapely fiona cartopy scikit-sparse suitesparse \
    pandas-plink statsmodels msprime xlrd openpyxl
conda activate feems
pip install git+https://github.com/NovembreLab/feems
```

Known problems (from the original analysis and from testing these scripts):

| Error | Fix |
|---|---|
| `ImportError: libpoppler.so...` when importing feems | `conda install -c conda-forge poppler` |
| `TypeError: 'tuple' object is not callable` in `_rank_one_solver` | scikit-sparse ≥ 0.5 is not compatible: `pip install "scikit-sparse<0.5"` |
| Slow or failing `conda install` | Create a new environment with only FEEMS; do not add it to an existing one |

## 3. Input files

- **VCF**: `run_feems.py` reads the genotypes directly from the VCF. No
  plink conversion is needed. Missing genotypes are replaced by the mean of
  the SNP, as in the FEEMS tutorials.
- **Coordinates**: `<sample ID> <tab> <longitude> <tab> <latitude>`, no
  header (as in [08](../08_construct/)).
- **Habitat polygon** (`.outer`): `longitude latitude`, one vertex per line,
  first point = last point. Make it as for EEMS
  ([09](../09_eems/#32-outer-the-habitat-polygon)), e.g. in QGIS or with
  <https://www.keene.edu/campus/maps/tool/>. A polygon with a few dozen
  vertices is enough; a polygon clipped to the coastline in QGIS can have
  tens of thousands of vertices, which slows everything down.

> **Longitude first.** In the original analysis, the first map was wrong
> because the coordinates were lat-long. FEEMS, like EEMS, needs
> **longitude latitude** in both files.

## 4. Build the grid

### 4.1 Why a custom grid?

FEEMS comes with three global grids:

| Grid | dggridR resolution | Distance between nodes |
|---|---|---|
| `grid_500` | 4 | ~440 km |
| `grid_250` | 5 | ~220 km |
| `grid_100` | 6 | ~110 km |

For a study area of a few hundred km, even `grid_100` places all samples on
a handful of nodes, and the surface says nothing. **This was the problem in
the original analysis** (populations within ~250 km): with `grid_100` and
`grid_250` the grid was too coarse, and FEEMS could not be used.

The solution is to make a finer grid for the study area only. The FEEMS
authors' rule of thumb: choose a resolution at which most individuals of
the same sampling site fall on one node, and different sites fall on
different nodes.

### 4.2 Make the grid

Edit the settings of `make_feems_grid.R`:

```r
OUTER     <- "feems/my_data.outer"
GRID_FILE <- "feems/grid_res8.shp"
RES       <- 8          # 6 ~110 km, 7 ~55 km, 8 ~28 km, 9 ~14 km
```

```bash
Rscript scripts/make_feems_grid.R
```

```
Resolution 8: cell spacing ~28 km, cell area ~389 km2
792 triangles written to feems/grid_res8.shp
```

Install dggridR with `install.packages("dggridR")`.

### 4.3 Check the graph

```bash
python3 scripts/run_feems.py --vcf filtered/my_data-miss9-maf1-th100.recode.vcf \
    --coords coords_individuals.txt --outer feems/my_data.outer \
    --grid feems/grid_res8.shp --outdir feems --prefix my_data --step graph
```

```
n_samples=51, n_snps=16064
Grid nodes: ..., observed nodes: ..., sampling sites: 11
```

Open `feems/my_data_graph.pdf`: the grid, the nodes with samples (observed
nodes) and the samples. If several sampling sites share a node, the script
prints a note; use a higher `RES`. Very fine grids (many nodes) are slower,
especially the cross-validation.

> **Maps.** The first time, cartopy downloads the coastlines from Natural
> Earth. On a server without internet, add `--no-map` (the plots then have
> no coastlines).

## 5. Cross-validation and fit

```bash
python3 scripts/run_feems.py --vcf filtered/my_data-miss9-maf1-th100.recode.vcf \
    --coords coords_individuals.txt --outer feems/my_data.outer \
    --grid feems/grid_res8.shp --outdir feems --prefix my_data --step all
```

```
Best lambda (lowest CV error): 0.04642  -> feems/my_data_cv.txt
Migration surface: feems/my_data_fit_lamb0.04642.pdf
```

(Numbers from the simulated test data.)

- **Cross-validation** (`_cv.pdf`, `_cv.txt`): FEEMS leaves out some
  observed nodes, fits the model on the rest and predicts the left-out
  nodes, for each λ (default 10 values from 0.001 to 100, 5 folds). The
  best λ has the lowest error. If it is at the edge of the range, widen it
  with `--lamb-min` / `--lamb-max`.
- **Fit** (`_fit_lamb<λ>.pdf`): the migration surface with the best λ. To
  see how sensitive the result is, fit with other values:

  ```bash
  python3 scripts/run_feems.py ... --step fit --lamb 1
  ```

## 6. Interpretation

- Orange edges (log weights below the mean) between groups of samples are
  candidate **barriers**; blue edges are **corridors**.
- Only interpret the surface between observed nodes. Edges far from the
  samples reflect the penalty, not the data.
- Compare with EEMS ([09](../09_eems/)): barriers found by both methods,
  and coinciding with clusters of [06](../06_pca_dapc/)–[10](../10_snmf_tess3/),
  are the most reliable.

## References

- Marcus J, Ha W, Barber RF, Novembre J (2021) Fast and flexible estimation
  of effective migration surfaces. *eLife* 10: e61927.
- FEEMS: <https://github.com/NovembreLab/feems> (see the notebooks
  `getting-started` and `cross-validation`).
- Barnes R, Sahr K (2017) dggridR: discrete global grids for R.
  <https://github.com/r-barnes/dggridR>
- Petkova D, Novembre J, Stephens M (2016) Visualizing spatial population
  structure with estimated effective migration surfaces. *Nature Genetics*
  48: 94–100.
