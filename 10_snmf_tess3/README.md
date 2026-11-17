# 10 — sNMF and TESS3: ancestry coefficients and maps

Estimate ancestry coefficients with two fast alternatives to STRUCTURE:
**sNMF** (LEA package, no geography) and **TESS3** (tess3r package, uses the
coordinates of the individuals), and map them as pie charts.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/snmf_tess3.R`](scripts/snmf_tess3.R) | Converts the VCF, runs sNMF and TESS3 for a range of K, plots to choose K, then Q tables, bar plots and pie maps |

**Requirements:** R with `vcfR`, `LEA`, `tess3r`, `ggplot2`, `maps` and
`scatterpie` (see below). A filtered VCF **without the outgroup** from
[04_filtering](../04_filtering/), a popmap and the coordinates of each
individual.

```r
BiocManager::install("LEA")
remotes::install_github("bcm-uga/TESS3_encho_sen")   # tess3r
install.packages(c("scatterpie", "maps"))
```

---

## 1. Background

Both methods estimate the same thing as STRUCTURE ([07](../07_structure/)):
for each individual, the proportion of its genome from each of K ancestral
populations (the Q matrix). They are much faster, so many values of K and
replicates can be run in minutes.

| | sNMF | TESS3 |
|---|---|---|
| Model | Matrix factorisation, no population genetic assumptions (e.g. no Hardy-Weinberg) | Same, plus a spatial term: nearby individuals tend to have similar ancestry |
| Uses coordinates | No | Yes |
| Choosing K | **Cross-entropy**: lowest = best | **Cross-validation score**: lowest = best |

sNMF does not assume Hardy-Weinberg equilibrium, so it is better than
STRUCTURE for inbred or selfing species. TESS3 is useful when structure is
spatial, but, like EEMS and conStruct ([08](../08_construct/),
[09](../09_eems/)), it can only be interpreted where there are samples.

## 2. Input files

- **VCF**, **popmap** (`<population> <tab> <sample ID>`).
- **Coordinates**: `<sample ID> <tab> <longitude> <tab> <latitude>`, no
  header (the same file as in [08](../08_construct/)).

> The coordinates are matched to the samples **by ID**. In the original
> analysis, the coordinates file had no IDs and was matched by order: if
> the order of the VCF changes (e.g. after removing a sample), all the
> coordinates shift without any error.

The script converts the VCF to the LEA `.geno` format itself (one line per
SNP, one digit per individual: 0, 1, 2 copies of the alternative allele,
9 = missing).

## 3. Run sNMF and TESS3 (pass 1)

Edit the settings of `snmf_tess3.R`:

```r
VCF      <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
POPMAP   <- "popmap.txt"
COORDS   <- "coords_individuals.txt"
K_MAX    <- 16
REPS     <- 10
ALPHA    <- 100
K_PLOT   <- NA
```

```bash
Rscript scripts/snmf_tess3.R
```

```
131 individuals, 13179 SNPs, ... populations
sNMF: K = 1..16, 10 repetitions. Start: ...
TESS3: K = 1..16, 10 repetitions. Start: ...
Choose K from snmf_tess3/*cross-entropy.pdf and *cross-validation.pdf, set K_PLOT and run again.
```

- `ALPHA` is the sNMF regularisation parameter. Higher values give smoother
  ancestry coefficients; the LEA documentation recommends increasing it
  for small datasets (fewer than ~100 individuals). 100 was used in the
  original analysis.
- The runs are saved (`.snmfProject`, `_tess3.rds`). The second pass loads
  them instead of running again (unless `RERUN <- TRUE`).

### Choose K

- `_snmf_cross-entropy.pdf` / `.txt`: cross-entropy for each K (mean and
  range of the replicates).
- `_tess3_cross-validation.pdf`: cross-validation score for each K.

Choose **K where the curve has its minimum or stops decreasing** (the
"elbow"). As with STRUCTURE, the curve often keeps decreasing slowly, so
look at the bar plots of a few values of K. In the original analysis
(131 individuals, K = 1–16), the cross-entropy suggested K = 7–9, and K = 9
was used for both sNMF and TESS3.

## 4. Tables and maps (pass 2)

```r
K_PLOT <- 9
```

```bash
Rscript scripts/snmf_tess3.R
```

For each method (sNMF uses the replicate with the lowest cross-entropy;
TESS3 the best run):

| File | Content |
|---|---|
| `_sNMF_K9_Q.csv`, `_TESS3_K9_Q.csv` | Q of each individual, its most likely cluster and the probability of that cluster, with population and coordinates |
| `_sNMF_K9_barplot.pdf` | Ancestry bar plot, individuals grouped by population |
| `_sNMF_K9_map.pdf` | Pie charts on the map: **one pie per site** (mean Q of the individuals at the same coordinates) |

The screen output also shows the table population × most likely cluster.

> **One pie per site.** With one pie per individual, the individuals of a
> site are drawn on top of each other and only the last one is visible. The
> script averages the individuals of each site. Adjust the size of the pies
> with `PIE_SCALE` and the map area with `MAP_MARGIN`.

**Colors.** The cluster numbers are arbitrary and differ between sNMF and
TESS3 (and between K). To use the same color for the same group in both
maps, reorder `COLORS`, or edit the final figure.

## 5. Interpretation

- **Compare the methods.** sNMF, TESS3, STRUCTURE and DAPC
  ([06](../06_pca_dapc/), [07](../07_structure/)) often agree on the main
  groups and differ in the details. Groups found by all methods are robust.
  In the original analysis, DAPC found 8–11 clusters and sNMF/TESS3 9.
- **TESS3 vs sNMF.** If TESS3 gives smoother transitions between nearby
  sites than sNMF, part of the structure is spatial (isolation by
  distance); test this with conStruct ([08](../08_construct/)).
- **Admixed individuals** are not necessarily hybrids (see
  [07](../07_structure/#6-interpretation)).

## References

- Frichot E, Mathieu F, Trouillon T, Bouchard G, François O (2014) Fast and
  efficient estimation of individual ancestry coefficients. *Genetics* 196:
  973–983. (sNMF)
- Frichot E, François O (2015) LEA: an R package for landscape and
  ecological association studies. *Methods in Ecology and Evolution* 6:
  925–929.
- Caye K, Deist TM, Martins H, Michel O, François O (2016) TESS3: fast
  inference of spatial population structure and genome scans for selection.
  *Molecular Ecology Resources* 16: 540–548.
- Tutorial used in the original analysis: Connor French, *Introduction to
  population structure in R*,
  <https://connor-french.github.io/intro-pop-structure-r/>
