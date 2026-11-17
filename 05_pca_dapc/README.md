# 05 — PCA and DAPC in RStudio

Explore the genetic structure without a population genetic model: a **PCA**
to see how the individuals are distributed, and a **DAPC** to find genetic
clusters and assign the individuals to them.

This tutorial is done **in RStudio, step by step**: some functions ask you
questions in the Console, and the answers depend on the plots you see.
Write your answers in your notes, so you can repeat the analysis.

## Install

In RStudio (once):

```r
install.packages(c("vcfR", "adegenet", "ggplot2"))
```

## Prepare the folder

```
05_pca_dapc/
├── data/        my_data-miss0.9-maf0.01-th100.recode.vcf  (from 03_filtering, without outgroup)
│                popmap.txt
└── results/     plots and tables
```

If the VCF is on a server, copy it to your computer (e.g. with FileZilla or
`scp`).

---

## 1. Background

- **PCA** (principal component analysis) summarises the genotypes in a few
  axes (principal components, PCs) that explain most of the variation. It
  does not use the population labels: the groups you see are in the data.
- **DAPC** (discriminant analysis of principal components) has two steps:
  1. `find.clusters` groups the individuals into K clusters (k-means on the
     PCs), without using the population labels;
  2. `dapc` builds discriminant functions that best **separate** these
     clusters, and gives each individual a probability of belonging to each
     cluster.

PCA shows the overall variation; DAPC maximises the differences between
groups. Always look at the PCA first.

> **Use the VCF without the outgroup.** An outgroup is so different from
> the other samples that it takes the first PC, and the structure within
> your species becomes hard to see.

## 2. Read the data

Open RStudio, create a new script (File → New File → R Script), and run each
block with Ctrl+Enter (Cmd+Enter on Mac).

```r
library(vcfR)
library(adegenet)
library(ggplot2)

setwd("~/my_project/05_pca_dapc")      # change to your folder

# read the VCF and convert it to a genlight object (adegenet format)
vcf <- read.vcfR("data/my_data-miss0.9-maf0.01-th100.recode.vcf")
gl  <- vcfR2genlight(vcf)

# add the population of each individual, from the popmap
popmap <- read.table("data/popmap.txt", header = FALSE, col.names = c("pop", "id"))
pop(gl) <- popmap$pop[match(indNames(gl), popmap$id)]

gl                 # summary: number of individuals, SNPs and populations
table(pop(gl))     # individuals per population
```

If `table(pop(gl))` shows `NA`, some sample IDs of the VCF are not in the
popmap (check the names, and Windows line endings, see
[00_setup](../00_setup/#21-files-used-in-many-tutorials)).

## 3. PCA

```r
pca <- glPca(gl)
```

R shows a bar plot of the eigenvalues (the variance explained by each PC)
and asks in the Console:

```
Select the number of axes:
```

Type the number of PCs to keep, e.g. `10`, and press Enter. The number only
affects how many PCs are saved (you will plot the first 2–3).

How much variance each PC explains:

```r
var_pc <- 100 * pca$eig / sum(pca$eig)
round(var_pc[1:5], 1)      # % of variance of PC1 to PC5
```

With many SNPs and closely related populations, PC1 often explains only a
few % of the variance: this is normal.

Plot PC1 × PC2, colored by population:

```r
pca_df <- data.frame(pca$scores, population = pop(gl))

ggplot(pca_df, aes(x = PC1, y = PC2, fill = population)) +
  geom_point(shape = 21, size = 3) +
  labs(x = paste0("PC1 (", round(var_pc[1], 1), "%)"),
       y = paste0("PC2 (", round(var_pc[2], 1), "%)")) +
  theme_bw()
```

Change `PC2` to `PC3` to see PC1 × PC3: the third axis can separate groups
that overlap on PC1–PC2 (in the original analysis, PC3 separated two
species that overlapped on PC1–PC2).

Save the plot with **Export → Save as Image / Save as PDF** in the Plots
panel, or with `ggsave("results/PCA_PC1_PC2.pdf")` after the plot.

> **Check:** an individual far from its population may be a mislabelled
> sample, a hybrid, or a sample with much missing data. Check its missing
> data ([03_filtering](../03_filtering/#4-filter)) before interpreting it.

The `dartR` package has similar functions (`gl.pcoa` and `gl.pcoa.plot`),
used in the original analysis.

## 4. DAPC step 1: find the clusters

```r
grp <- find.clusters(gl, max.n.clust = 15)
```

`max.n.clust` is the largest number of clusters tested. Set it a bit higher
than the number of groups you expect.

### Question 1: number of PCs

R shows a plot of the **cumulative variance** explained by the PCs and asks:

```
Choose the number PCs to retain (>= 1):
```

At this step, **keep all the PCs**: type a number equal to or larger than
the number of individuals (e.g. `150` for 135 individuals). The DAPC
tutorial also keeps all PCs at this step: no information is lost.

### Question 2: number of clusters (K)

R shows the **BIC** for each K and asks:

```
Choose the number of clusters (>=2):
```

How to choose:

- The best K has the **lowest BIC**. Often the BIC decreases and then
  increases again: choose the K at the minimum.
- If the BIC keeps decreasing slowly, choose the K at the **elbow** of the
  curve, where it stops decreasing strongly. A large K splits populations
  into small clusters with no biological meaning.
- If there is no clear elbow, try two or three values of K and compare the
  results. As the DAPC tutorial says, the question is "how many clusters
  are useful to describe the data?", not the "true" K.

Then see how the populations fall into the clusters:

```r
table(pop(gl), grp$grp)
table.value(table(pop(gl), grp$grp), col.lab = paste("cluster", 1:length(levels(grp$grp))))
```

Rows are populations, columns are clusters: does each population fall into
one cluster? Are some merged, or split into several clusters?

## 5. DAPC step 2: describe the clusters

### 5.1 How many PCs? (a-score)

The DAPC uses some PCs to build the discriminant functions. **With too many
PCs, the DAPC over-fits**: it separates any groups perfectly, even random
ones, and the membership probabilities become unreliable. The **a-score**
helps to choose: it compares the success of the assignment with the success
for random groups.

Run a first DAPC with many PCs, answering the questions with all PCs and
all discriminant functions:

```r
dapc_test <- dapc(gl, grp$grp)
```

```
Choose the number PCs to retain (>=1):          <- type the number of individuals, e.g. 135
Choose the number discriminant functions to retain (>=1):   <- type K - 1
```

Then:

```r
a <- optim.a.score(dapc_test)
a$best
```

The plot shows the a-score for different numbers of PCs; the red line is
the best value (`a$best`). In the original analysis: K = 8 and 12 PCs.

### 5.2 Final DAPC

Run the DAPC again:

```r
dapc1 <- dapc(gl, grp$grp)
```

```
Choose the number PCs to retain (>=1):          <- the a-score value, e.g. 12
Choose the number discriminant functions to retain (>=1):   <- K - 1 (all), e.g. 7
```

With few clusters, keep all the discriminant functions (K − 1). With many
clusters (tens), keep the first few.

To repeat the analysis without the questions, give the answers as
arguments:

```r
grp   <- find.clusters(gl, max.n.clust = 15, n.pca = 150, n.clust = 8)
dapc1 <- dapc(gl, grp$grp, n.pca = 12, n.da = 7)
```

`find.clusters` starts from random positions, so the cluster numbers can
change between runs. Run `set.seed(123)` before it to always get the same
result.

## 6. DAPC plots

**Scatter plot** of the individuals on the first two discriminant
functions, with the clusters as ellipses:

```r
scatter(dapc1, posi.da = "bottomright", bg = "white", pch = 20, cstar = 0,
        col = funky(8), scree.pca = TRUE, posi.pca = "bottomleft",
        leg = TRUE, txt.leg = paste("Cluster", 1:8))
```

The small inset plots show the eigenvalues of the discriminant functions and
the PCs used.

**Membership probabilities** (a STRUCTURE-like bar plot):

```r
compoplot(dapc1, posi = "bottomright", txt.leg = paste("Cluster", 1:8),
          lab = "", ncol = 1, col = funky(8), show.lab = FALSE)
```

**Table** with the cluster and membership probabilities of each individual:

```r
memb <- data.frame(individual = indNames(gl), population = pop(gl),
                   cluster = dapc1$assign, round(dapc1$posterior, 3))
write.csv(memb, "results/DAPC_K8_membership.csv", row.names = FALSE)
```

## 7. Interpretation

- **Population × cluster table.** Compare the clusters with the taxonomy
  and the geography.
- **Membership probabilities** are usually close to 0 or 1. An individual
  with intermediate probabilities (e.g. 0.6/0.4) is genetically
  intermediate, but **DAPC is not a test for hybrids**: use STRUCTURE
  ([06](../06_structure/)) and the introgression tutorials
  ([phylogenomics](../phylogenomics/)).
- **Robustness.** Run the analysis with VCFs filtered for different amounts
  of missing data ([03](../03_filtering/)). Clusters that remain are
  robust.

## Next step

Model-based clustering: [06_structure](../06_structure/).

## More information

- **DAPC tutorial** by Thibaut Jombart (detailed explanation of each step,
  the a-score, cross-validation with `xvalDapc`, and how to customise the
  plots): <https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf>
- adegenet tutorials: <https://adegenet.r-forge.r-project.org/>
- Introduction to population structure in R (DAPC and sNMF, with maps):
  <https://connor-french.github.io/intro-pop-structure-r/>

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
