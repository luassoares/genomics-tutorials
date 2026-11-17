# 08 — conStruct: structure and isolation by distance

Test whether genetic differentiation is better explained by discrete
clusters, by continuous isolation by distance, or by both, with conStruct.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/run_construct.R`](scripts/run_construct.R) | Prepares the data from a VCF and coordinates, runs the cross-validation and the full spatial and non-spatial runs |
| [`scripts/plot_construct.R`](scripts/plot_construct.R) | Cross-validation plots and t-tests, layer contributions, admixture bar plot and pie map |

**Requirements:** R with `conStruct`, `vcfR`, `geodist`, `doParallel` and
`foreach`. A filtered VCF **without the outgroup** from
[04_filtering](../04_filtering/), and the coordinates of each individual.

---

## 1. Background

STRUCTURE-like methods ([07](../07_structure/)) assume discrete
populations. When populations are distributed along a continuous range,
**isolation by distance** (IBD: nearby individuals are more similar) can
appear as several clusters with "admixed" individuals between them, even
without discrete groups.

conStruct models each individual as a mix of **K layers**. In the
**spatial model**, relatedness within each layer decreases with geographic
distance; the **non-spatial model** is like STRUCTURE. Comparing them
answers two questions:

1. Does adding distance improve the model? (spatial vs non-spatial)
2. How many layers are needed? (cross-validation and layer contributions)

## 2. Input files

- **VCF**: the script calculates the allele frequency of each individual
  (0, 0.5, 1) directly from the genotypes.
- **Coordinates**: one line per individual, `<sample ID> <tab> <longitude>
  <tab> <latitude>`, decimal degrees, **longitude first**:

  ```
  1023	-50.950	-25.200
  1024	-50.950	-25.200
  1101	-51.217	-25.333
  ```

  Individuals from the same site can have the same coordinates.

> **Missing data.** Individuals with much missing data make the MCMC
> unstable. Remove them with `EXCLUDE` or `MAX_MISSING_IND`. In the
> original analysis, five individuals with much missing data were removed
> before the runs mixed well.

## 3. Run

Edit the settings of `run_construct.R`:

```r
VCF       <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
COORDS    <- "coords_individuals.txt"
K_MAX     <- 5
XVAL_REPS <- 10
XVAL_ITER <- 1e4
N_CORES   <- 16
```

```bash
Rscript scripts/run_construct.R
```

```
80 samples, 8120 loci
Cross-validation: K = 1..5, 10 replicates. Start: ...
```

The script does two steps (each can be switched off with `RUN_XVAL` /
`RUN_FULL`):

1. **Cross-validation** (`x.validation`): for each replicate, 90% of the
   loci are used to fit the model and the other 10% to measure how well it
   predicts them (predictive accuracy), for spatial and non-spatial
   models, K = 1–5. This is the slow part: **several hours** with 16 cores
   for ~80 individuals.
2. **Full runs** of each model and K, for the layer contributions and the
   plots.

> **Check the MCMC.** conStruct uses Stan. Warnings such as `divergent
> transitions`, `R-hat > 1.05` or `Effective Samples Size (ESS) is too low`
> mean the chain has not converged: increase `FULL_ITER`, and check for
> individuals with much missing data. Look at the trace plots written by
> each run (`*trace.plots*.pdf`).

## 4. Results

```r
OUTDIR <- "construct"
K_MAX  <- 5
K_PLOT <- 3
```

```bash
Rscript scripts/plot_construct.R
```

### 4.1 Spatial or non-spatial? How many layers?

`_xval.pdf` shows the predictive accuracy (higher = better) for each K and
model. `_xval_ttests.txt` has paired t-tests:

```
Spatial K=2 vs K=1: mean difference = 901.92, t = 33.20, p = 5.03e-11
...
K=1 spatial vs non-spatial: mean difference = ..., p = ...
```

- **Spatial > non-spatial** (especially at low K): geographic distance
  explains part of the differentiation (IBD). In the original analysis the
  spatial models had clearly higher predictive accuracy.
- **Choosing K**: the predictive accuracy usually increases quickly and
  then levels off. Choose the smallest K after which adding layers gives
  little improvement. A significant t-test is not enough: with many loci,
  tiny improvements are significant.

> The conStruct vignette computes the confidence interval with
> `sd(x)/length(x)`. The standard error is `sd(x)/sqrt(length(x))`, used in
> `plot_construct.R`; the intervals are therefore wider than in the
> vignette.

### 4.2 Layer contributions

`_layer_contributions.pdf`: for each K, how much each layer contributes
to the total covariance. A layer that contributes very little (e.g. < 1%)
is not needed: the best K is the number of layers with a real
contribution.

### 4.3 Admixture plots

For `K_PLOT`:

- `_spK3_structure.pdf`: proportions of each layer per individual, like a
  STRUCTURE plot.
- `_spK3_pies.pdf`: the same proportions as pie charts on the map.

## 5. Interpretation

- If the spatial model with **K = 1** is about as good as higher K, the
  structure is mostly continuous (IBD), and the clusters of
  STRUCTURE/DAPC may be artefacts of IBD.
- If several spatial layers are needed, there are discrete groups **in
  addition to** IBD.
- Compare the layers with the clusters of [06](../06_pca_dapc/) and
  [07](../07_structure/), and with the migration barriers of EEMS
  ([09](../09_eems/)).

## References

- Bradburd GS, Coop GM, Ralph PL (2018) Inferring continuous and discrete
  population genetic structure across space. *Genetics* 210: 33–52.
- conStruct vignettes: `vignette(package = "conStruct")` — format-data,
  run-conStruct, visualize-results, model-comparison.
