# 09 — EEMS: effective migration surfaces

Map where gene flow is lower (barriers) or higher (corridors) than expected
under isolation by distance, with EEMS.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/vcf2diffs.R`](scripts/vcf2diffs.R) | Makes the `.diffs` and `.coord` input files from a VCF (no plink or bed2diffs needed) |
| [`scripts/run_eems.sh`](scripts/run_eems.sh) | Writes the parameter files and runs several EEMS chains in parallel |
| [`scripts/plot_eems.R`](scripts/plot_eems.R) | Migration and diversity maps, model-fit plots (reemsplots2) |

**Requirements:** EEMS (`runeems_snps`, compiled below); R with `vcfR`,
`reemsplots2`, `ggplot2` and `maps`. A filtered VCF **without the outgroup**
from [04_filtering](../04_filtering/), and the coordinates of each
individual.

---

## 1. Background

EEMS places a grid of **demes** over the study area and estimates two
surfaces:

- **Effective migration (m)**: how fast genetic similarity decays between
  nearby locations. Areas with **lower** m than average are **barriers**
  to gene flow (orange in the plots); areas with higher m are corridors
  (blue).
- **Effective diversity (q)**: the expected genetic dissimilarity between
  two individuals from the **same** location (within-deme diversity).

Without barriers or corridors (pure isolation by distance), the migration
surface is flat. Colors are on a log10 scale relative to the mean: a value
of 1 means 10 times more migration than average, −1 means 10 times less.

## 2. Install EEMS

```bash
cd ~/programs
git clone https://github.com/dipetkov/eems.git
cd eems/runeems_snps/src
```

EEMS needs the **Eigen** and **Boost** libraries. Eigen is header-only:

```bash
cd ~/programs
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xzf eigen-3.4.0.tar.gz
```

Boost is often already installed on servers (look for `/usr/include/boost`).
If not, install it with conda: `conda install -c conda-forge boost`.

Edit the three paths at the top of the `Makefile`, then compile:

```bash
cd ~/programs/eems/runeems_snps/src
nano Makefile
```

```
EIGEN_INC = /home/<user>/programs/eigen-3.4.0
BOOST_LIB = /usr/lib/x86_64-linux-gnu        # or $CONDA_PREFIX/lib
BOOST_INC = /usr/include                     # the folder that contains boost/ (or $CONDA_PREFIX/include)
```

```bash
make linux
./runeems_snps --help
echo 'export PATH=$HOME/programs/eems/runeems_snps/src:$PATH' >> ~/.bashrc
```

> **You do not need bed2diffs.** The EEMS README converts the data with
> `bed2diffs`, which needs plink files and the `libplinkio` library (often
> hard to compile). `vcf2diffs.R` computes the same matrix directly from
> the VCF, with the formulas of the EEMS `str2diffs` documentation.

## 3. Input files

EEMS needs three files with the same prefix (`datapath`):

| File | Content |
|---|---|
| `<datapath>.diffs` | Matrix of average genetic differences between all pairs of individuals |
| `<datapath>.coord` | Coordinates of each individual, **longitude latitude**, same order as `.diffs` |
| `<datapath>.outer` | Polygon of the habitat, **longitude latitude**, counterclockwise, first point = last point |

### 3.1 `.diffs` and `.coord`

Coordinates file: `<sample ID> <tab> <longitude> <tab> <latitude>`, decimal
degrees (the same file as in [08_construct](../08_construct/)).

```r
VCF      <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
COORDS   <- "coords_individuals.txt"
DATAPATH <- "eems/data/my_data"
```

```bash
Rscript scripts/vcf2diffs.R
```

```
Written eems/data/my_data.diffs and eems/data/my_data.coord
nIndiv = 131
nSites = 13179
Now make eems/data/my_data.outer (habitat polygon, see the tutorial).
```

Note `nIndiv` and `nSites` for the parameter file.

`METHOD <- "v1"` uses, for each pair of individuals, only the SNPs
genotyped in both. `"v2"` replaces missing genotypes by the mean. If the
script warns that the matrix has more than one positive eigenvalue (not a
valid dissimilarity matrix), use `"v2"`.

### 3.2 `.outer`: the habitat polygon

The polygon outlines the area where the species can occur. Draw it around
the sampling points with some margin, and follow natural limits (e.g. the
coast), so that EEMS does not place demes in the sea.

- Online: draw the polygon on a map with a tool such as
  <https://www.birdtheme.org/useful/v3tool.html> and copy the coordinates.
- QGIS: draw a polygon, clip it with a country or coastline shapefile, and
  export the vertices.

The file has one vertex per line, **longitude latitude**, listed
**counterclockwise**, and the last line repeats the first:

```
-57.97 -20.40
-57.79 -32.06
-44.30 -31.95
-44.13 -20.20
-57.97 -20.40
```

> **Longitude first, in all three files.** Mixing the order (lat-long in
> one file, long-lat in another) is the most common error: the samples fall
> outside the polygon or the map appears rotated.

## 4. Run EEMS

Edit the settings of `run_eems.sh`:

```bash
DATAPATH=eems/data/my_data
NINDIV=131
NSITES=13179
NDEMES=300
NCHAINS=3
SEEDS="123 456 789"
NUM_MCMC=4000000
NUM_BURN=1000000
NUM_THIN=9999
```

```bash
bash scripts/run_eems.sh
```

- **nDemes**: density of the grid. Few demes (e.g. 100–200) run faster; more
  demes (300–500) give finer maps. Try two densities: real barriers appear
  in both.
- **Several chains** (different seeds) check convergence: the maps of all
  chains should be similar. `plot_eems.R` averages them.

## 5. Plot

```r
MCMCPATH <- Sys.glob("eems/runs/my_data-nDemes300-chain*")
COORD    <- "eems/data/my_data.coord"
```

```bash
Rscript scripts/plot_eems.R
```

Install reemsplots2 once with `remotes::install_github("dipetkov/reemsplots2")`.

### 5.1 Check the runs first

- **`_pilogl01.pdf`**: log posterior of each chain along the MCMC. The
  chains should reach the same level and stay flat. If they are still
  increasing, run longer (`NUM_MCMC`, `NUM_BURN`).
- **`_rdist01.pdf` and `_rdist02.pdf`**: observed vs fitted dissimilarities
  between and within demes. A **strong linear relationship** means the model
  fits the data well.
- **`_rdist03.pdf`**: genetic dissimilarity vs geographic distance. Under
  isolation by distance it increases linearly with distance.

### 5.2 Migration and diversity maps

- **`_mrates01.pdf`** (effective migration): orange = lower migration than
  average (barrier); blue = higher (corridor). Circles = sampling sites.
- **`_mrates02.pdf`**: posterior probability that migration is lower or
  higher than average. Only interpret areas where it is high (e.g. > 0.9).
- **`_qrates01.pdf` / `_qrates02.pdf`**: effective diversity. In the
  original analysis, inland populations had higher diversity than
  populations at the edges of the range.

## 6. Interpretation and pitfalls

- **Only interpret areas with samples.** Far from the sampling sites the
  surface reflects the prior, not the data.
- **Barriers between populations that are also separated in
  STRUCTURE/DAPC** ([06](../06_pca_dapc/), [07](../07_structure/)) support
  real limits to gene flow. Compare them with landscape features
  (rivers, valleys, changes in altitude or habitat).
- **FEEMS** is a faster successor of EEMS. Its default grids are too
  coarse for small study areas, but a finer grid can be made: see
  [11_feems](../11_feems/).

## References

- Petkova D, Novembre J, Stephens M (2016) Visualizing spatial population
  structure with estimated effective migration surfaces. *Nature Genetics*
  48: 94–100.
- EEMS: <https://github.com/dipetkov/eems>; reemsplots2:
  <https://github.com/dipetkov/reemsplots2>
- Marcus J, Ha W, Barber RF, Novembre J (2021) Fast and flexible estimation
  of effective migration surfaces. *eLife* 10: e61927. (FEEMS)
