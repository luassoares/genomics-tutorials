# 07 — STRUCTURE and fastSTRUCTURE

Estimate the number of genetic clusters (K) and the ancestry proportions of
each individual, with STRUCTURE and fastSTRUCTURE run through
structure_threader.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/vcf2structure.py`](scripts/vcf2structure.py) | Converts a VCF to STRUCTURE and fastSTRUCTURE input, ordered by population, plus the population and label files |
| [`scripts/plot_structure.R`](scripts/plot_structure.R) | Evanno table, alignment of the replicates, bar plots with population labels (pophelper) |

**Requirements:** STRUCTURE, fastSTRUCTURE and structure_threader; R with
`pophelper`; Python 3. A filtered VCF **without the outgroup** from
[04_filtering](../04_filtering/).

Install structure_threader with `pip install structure_threader --user`. It
includes STRUCTURE and fastSTRUCTURE binaries; their paths are given with
`-st` and `-fs` (see the [structure_threader
documentation](https://structure-threader.readthedocs.io/)).

---

## 1. Background

STRUCTURE assumes that there are K ancestral populations (clusters), each
with its own allele frequencies, and estimates for each individual the
proportion of its genome that comes from each cluster (Q). It is a Bayesian
model run with MCMC, so:

- each K is run several times (**replicates**) to check that the runs agree;
- K is not estimated directly: you compare the runs for different K.

**fastSTRUCTURE** uses the same model with a much faster algorithm
(variational inference). It is useful to explore many K quickly; STRUCTURE
is the reference method. Running both and comparing them is a good check.

## 2. Input files

```bash
python3 scripts/vcf2structure.py \
    --vcf filtered/my_data-miss9-maf1-th100.recode.vcf \
    --popmap popmap.txt --out my_data
```

```
Individuals (NUMINDS): 135
Loci (NUMLOCI):        13179
Populations:           13
Written: my_data.str, my_data_fast.str, my_data_popfile.txt, my_data_individuals.txt
```

| File | Used by |
|---|---|
| `my_data.str` | STRUCTURE |
| `my_data_fast.str` | fastSTRUCTURE |
| `my_data_popfile.txt` | structure_threader `--pop` (population names in its plots) |
| `my_data_individuals.txt` | `plot_structure.R` (population labels) |

The individuals are ordered by population (in the order of the popmap), so
the bar plots show each population as a block.

> **Other converters.** PGDSpider or `dartR::gl2structure` also work, but
> the columns of their output differ: check that the settings in
> `mainparams` (`LABEL`, `POPDATA`, `EXTRACOLS`, `MARKERNAMES`) match the
> file. If a file was made on Windows, convert it first with `dos2unix`.

## 3. STRUCTURE

### 3.1 Parameter files

```bash
mkdir structure && cd structure
structure_threader params -o .        # writes mainparams and extraparams
```

Edit `mainparams`. The important lines, for the files made by
`vcf2structure.py`:

```
#define BURNIN    125000   // (int) length of burnin period
#define NUMREPS   500000   // (int) number of MCMC reps after burnin

#define NUMINDS    135     // number of individuals (printed by vcf2structure.py)
#define NUMLOCI    13179   // number of loci (printed by vcf2structure.py)
#define PLOIDY       2
#define MISSING     -9
#define ONEROWPERIND 0     // two rows per individual

#define LABEL        1     // first column: individual ID
#define POPDATA      1     // second column: population number
#define POPFLAG      0
#define LOCDATA      0
#define PHENOTYPE    0
#define EXTRACOLS    0
#define MARKERNAMES  1     // first row: marker names
```

`MAXPOPS`, `INFILE` and `OUTFILE` are set by structure_threader. Leave
`extraparams` with the defaults (admixture model, correlated allele
frequencies). `POPDATA 1` only reads the population column: it is not used
in the analysis unless `USEPOPINFO` is set to 1 in `extraparams`.

> **Burn-in and run length.** 125,000 burn-in and 500,000 MCMC steps are
> safe values for SNP data. Shorter runs (e.g. 10,000 / 50,000) are fine to
> test the pipeline. Check convergence: the replicates of the same K should
> have similar likelihoods.

### 3.2 Run

```bash
structure_threader run -K 15 -R 10 -i ../my_data.str -o structure_out \
    -t 25 -st ~/.local/bin/structure --pop ../my_data_popfile.txt
```

- `-K 15`: run K = 1 to 15.
- `-R 10`: 10 replicates per K.
- `-t 25`: runs in parallel (each run uses one CPU).
- `-st`: path of the STRUCTURE binary.

15 K × 10 replicates = 150 runs. With 135 individuals and 13,000 SNPs this
can take days: run it inside `byobu`/`tmux` (see
[00_setup](../00_setup/#4-long-runs)).

## 4. fastSTRUCTURE

```bash
structure_threader run -K 15 -R 10 -i ../my_data_fast.str -o faststructure_out \
    -t 15 -fs ~/.local/bin/fastStructure --pop ../my_data_popfile.txt
```

It runs in minutes to hours. structure_threader also runs fastSTRUCTURE's
`chooseK`:

```bash
cat faststructure_out/chooseK.txt
```

```
Model complexity that maximizes marginal likelihood = 2
Model components used to explain structure in data = 2
```

The two values give a range of plausible K.

## 5. Choose K and plot

Edit the settings of `plot_structure.R`:

```r
RESULTS_DIR <- "structure/structure_out"
PROGRAM     <- "structure"            # or "faststructure"
INDIVIDUALS <- "my_data_individuals.txt"
K_PLOT      <- 2:6
```

```bash
Rscript scripts/plot_structure.R
```

For STRUCTURE the script writes the **Evanno table** (`_evanno.txt`) and
plots (`evannoMethodStructure.pdf`): mean ln P(D) for each K, and ΔK.
Then it aligns the replicates of each K (the cluster numbers can differ
between runs), averages them, and plots K = 2–6 in one figure, with
population labels.

**How to choose K:**

- **ΔK** (Evanno) has a peak at the K where the likelihood increases most.
  It often gives **K = 2**, because it detects the uppermost level of
  structure (e.g. two major lineages). This does not mean that there is no
  further structure.
- **ln P(D)** usually increases and then reaches a plateau; the start of
  the plateau is another estimate of K.
- **Look at the bar plots for several K.** The "best" K is a model choice,
  not a biological truth. Report the K suggested by the methods, and
  describe what higher K values show (e.g. at K = 3 a third species
  separates).
- **Uneven sampling** (some populations with many more individuals) biases
  ΔK. Puechmaille (2016) proposes alternative estimators.

> **Example from the original analysis.** STRUCTURE and fastSTRUCTURE both
> chose K = 2, splitting the samples into two groups of species. At K = 3
> one more species separated. A DAPC ([06](../06_pca_dapc/)) found more
> clusters. The methods are not contradictory: they describe structure at
> different levels.

## 6. Interpretation

- **Admixed individuals** (bars with two or more colors) can be hybrids,
  but mixed ancestry can also come from incomplete lineage sorting or from
  a population that is not in the sample. Confirm with ancestry painting
  ([02](../02_ancestry_painting/)) and D statistics ([03](../03_dsuite/)).
- **Isolation by distance** produces gradients of ancestry along geography,
  which STRUCTURE can show as admixture. If populations are spread along a
  gradient, run conStruct ([08](../08_construct/)).
- **Missing data:** individuals with much missing data can have
  intermediate Q. Check them (see [04](../04_filtering/)).

## References

- Pritchard JK, Stephens M, Donnelly P (2000) Inference of population
  structure using multilocus genotype data. *Genetics* 155: 945–959.
- Falush D, Stephens M, Pritchard JK (2003) Inference of population
  structure using multilocus genotype data: linked loci and correlated
  allele frequencies. *Genetics* 164: 1567–1587.
- Raj A, Stephens M, Pritchard JK (2014) fastSTRUCTURE: variational
  inference of population structure in large SNP data sets. *Genetics* 197:
  573–589.
- Pina-Martins F, Silva DN, Fino J, Paulo OS (2017) Structure_threader: an
  improved method for automation and parallelization of programs
  structure, fastStructure and MavericK on multicore CPU systems. *Molecular
  Ecology Resources* 17: e268–e274.
- Evanno G, Regnaut S, Goudet J (2005) Detecting the number of clusters of
  individuals using the software STRUCTURE: a simulation study. *Molecular
  Ecology* 14: 2611–2620.
- Francis RM (2017) pophelper: an R package and web app to analyse and
  visualize population structure. *Molecular Ecology Resources* 17: 27–32.
- Puechmaille SJ (2016) The program STRUCTURE does not reliably recover the
  correct population structure when sampling is uneven: subsampling and new
  estimators alleviate the problem. *Molecular Ecology Resources* 16:
  608–627.
