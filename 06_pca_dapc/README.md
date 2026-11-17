# 06 — PCA and DAPC

Explore genetic structure without a model: a PCA to see how individuals are
distributed, and a DAPC to find genetic clusters and assign individuals to
them.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/pca_dapc.R`](scripts/pca_dapc.R) | PCA plots, choice of K (BIC), choice of the number of PCs (a-score), DAPC plots and membership table |

**Requirements:** R with the packages `vcfR`, `adegenet` and `ggplot2`. A
filtered VCF **without the outgroup** from [04_filtering](../04_filtering/).

---

## 1. Background

- **PCA** summarises the genotypes in a few axes (principal components)
  that explain most of the variation. It makes no assumptions about
  populations: the groups you see are in the data.
- **DAPC** (discriminant analysis of principal components) first finds
  clusters of individuals with k-means on the PCs (`find.clusters`), then
  builds discriminant functions that best **separate** those clusters. It
  gives each individual a probability of membership in each cluster.

PCA shows the overall variation; DAPC maximises the differences between
groups. Always look at the PCA first.

> **Use the VCF without the outgroup.** An outgroup is so different from the
> ingroup that it takes the first PC, and the structure within the ingroup
> becomes hard to see.

## 2. Settings

```r
VCF        <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
POPMAP     <- "popmap.txt"
PREFIX     <- "my_data"
MAX_K      <- 15
K          <- NA
N_PCA_DAPC <- NA
```

The script is run three times: each run stops after making the plot you
need to choose the next setting. This replaces the interactive questions of
`find.clusters` ("Choose the number PCs to retain", "Choose the number of
clusters"), so the analysis can be repeated with exactly the same settings.

## 3. PCA (pass 1)

```bash
Rscript scripts/pca_dapc.R
```

```
135 individuals, 13179 SNPs, 13 populations
PCA plots written to pca_dapc/
```

Look at:

- `_PCA_variance.svg`: how much variation each PC explains. With many SNPs
  and closely related populations, PC1 often explains only a few %: this is
  normal.
- `_PCA_PC1-PC2.svg` and `_PCA_PC1-PC3.svg`: individuals colored by
  population. PC3 can separate groups that overlap on PC1–PC2 (in the
  example data, PC3 separated two species that overlapped on PC1–PC2).

> **Check:** an individual far away from its population may be a
> mislabelled sample, a hybrid, or a sample with much missing data. Check
> its missing data (`vcftools --missing-indv`, [04](../04_filtering/))
> before interpreting it.

## 4. DAPC

### 4.1 Choose K (pass 1)

The same run also writes `_DAPC_BIC.svg` and prints the BIC for each K:

```
 K  BIC
 1  ...
 2  ...
 ...
```

Choose **K at the elbow of the curve**: where the BIC stops decreasing
strongly. The lowest BIC is not always the best choice: the curve often
keeps decreasing slowly, and a large K splits populations into small,
meaningless clusters. If there is no clear elbow, compare two or three
values of K.

`find.clusters` uses all PCs by default (`N_PCA_CLUST <- NA`), as the
adegenet authors recommend: at this step, using more PCs does not cause
over-fitting.

### 4.2 Choose the number of PCs (pass 2)

Set K and run again:

```r
K <- 8
```

```bash
Rscript scripts/pca_dapc.R
```

```
a-score: best number of PCs = 12. Set N_PCA_DAPC (e.g. 12) and run again.
```

This step is important: with **too many PCs, the DAPC over-fits** and
separates any groups perfectly, even random ones. The a-score
(`optim.a.score`) measures how much better the assignment is than with
random groups, and suggests the number of PCs that maximises it (red line
in `_DAPC_ascore.svg`).

### 4.3 Final DAPC (pass 3)

```r
N_PCA_DAPC <- 12
```

```bash
Rscript scripts/pca_dapc.R
```

(Illustrative output: population C is split between clusters 4 and 8.)

```
Populations (rows) x DAPC clusters (columns):
          cluster
population  1  2  3  4  5  6  7  8
     PopA  12  0  0  0  0  0  0  0
     PopB   0 23  0  0  0  0  0  0
     PopC   0  0  0  7  0  0  0  3
...
```

Output:

| File | Content |
|---|---|
| `_DAPC_K8_scatter.svg` | Individuals on the first two discriminant functions, with the eigenvalues of the PCA and DA |
| `_DAPC_K8_compoplot.svg` | Membership probabilities of each individual, ordered by population |
| `_DAPC_K8_membership.csv` | Cluster and membership probabilities of each individual |

## 5. Interpretation

- **Population × cluster table.** Does each population fall into one
  cluster? Are some populations merged into one cluster, or split into
  several? Compare it with the taxonomy and geography.
- **Membership probabilities.** In DAPC they are usually close to 0 or 1.
  An individual with intermediate probabilities (e.g. 0.6/0.4) is
  genetically intermediate, but **DAPC is not a test for hybrids**: use
  STRUCTURE ([07](../07_structure/)), ancestry painting
  ([02](../02_ancestry_painting/)) or D statistics ([03](../03_dsuite/)).
- **Robustness.** Run the analysis with VCFs filtered for different amounts
  of missing data. Clusters that remain are robust.

## References

- Jombart T, Devillard S, Balloux F (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. *BMC Genetics* 11: 94.
- Jombart T, Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of
  genome-wide SNP data. *Bioinformatics* 27: 3070–3071.
- Knaus BJ, Grünwald NJ (2017) vcfR: a package to manipulate and visualize
  variant call format data in R. *Molecular Ecology Resources* 17: 44–53.
- Miller JM, Cullingham CI, Peery RM (2020) The influence of a priori
  grouping on inference of genetic clusters: simulation study and
  literature review of the DAPC method. *Heredity* 125: 269–280.
