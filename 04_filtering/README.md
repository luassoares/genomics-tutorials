# 04 — Filtering the VCF and removing outlier SNPs

Prepare the VCFs used in the other tutorials: select the samples and keep
good-quality biallelic SNPs.

**Scripts**

| Script | What it does |
|---|---|
| [`scripts/filter_vcf.sh`](scripts/filter_vcf.sh) | Selects samples, removes monomorphic sites, filters with several missing-data thresholds and counts the SNPs |
| [`scripts/pcadapt_outliers.R`](scripts/pcadapt_outliers.R) | Optional, not recommended as a default filter: finds outlier SNPs with pcadapt (section 4) |

**Requirements:** bcftools, vcftools, R with the packages `pcadapt` and
`vcfR` (see [00_setup](../00_setup/)).

---

## 1. Which VCF for which analysis?

Each analysis needs slightly different data, so you will usually make more
than one VCF from the same raw file:

| VCF | Filters | Used in |
|---|---|---|
| **Thinned** | biallelic, missing, MAF, `--thin` | PCA, DAPC, STRUCTURE, conStruct, EEMS, SVDQuartets, diversity |
| **Not thinned** | biallelic, missing, low MAF | Dsuite ([03](../03_dsuite/)), ancestry painting ([02](../02_ancestry_painting/)) |
| **With outgroup** | as above, with the outgroup samples | SVDQuartets, Dsuite (needs an outgroup) |
| **Without outgroup** | as above, without the outgroup | PCA, DAPC, STRUCTURE (the outgroup would dominate the first axes) |

Why these filters:

- **Biallelic SNPs only** (`--min-alleles 2 --max-alleles 2 --remove-indels`):
  most methods assume two alleles per site.
- **Missing data** (`--max-missing 0.9` = a site must be genotyped in at
  least 90% of the samples): too much missing data adds noise and can create
  false clusters of individuals with many missing genotypes.
- **MAF** (`--maf 0.01`): removes very rare alleles, which are often
  sequencing errors and carry little information about structure.
- **Thinning** (`--thin 100` = at most one SNP per 100 bp): with RAD/DArT
  loci, SNPs on the same read are linked. Most structure and tree methods
  assume independent SNPs, so keep one SNP per locus.

## 2. Select the samples

Start from the VCF with all samples and keep only the ones for this project.
Make a list with one sample ID per line, exactly as in the VCF:

```bash
bcftools query -l all_samples.vcf | head      # see the sample IDs in the VCF
```

```
1023
1024
1101
```

**Replicates and duplicates.** If an individual was sequenced twice, keep
only the replicate with less missing data. Compare them with:

```bash
vcftools --vcf all_samples.vcf --missing-indv --out all_samples
sort -k5,5g all_samples.imiss | head        # 5th column: proportion of missing genotypes
```

> **Tip:** after removing samples, some sites have only one allele left
> (monomorphic). Remove them with
> `bcftools view -e 'AC==0 || AC==AN'`. `filter_vcf.sh` does this
> automatically. If you remove samples later (as when a duplicate is found
> after filtering), run this command again.

## 3. Filter

Edit the settings at the top of `filter_vcf.sh`:

```bash
INPUT=all_samples.vcf
PREFIX=my_data
SAMPLES=list_my_samples.txt
EXCLUDE=""                       # e.g. "3419833,3419842" to drop duplicates
MISSING="0.5 0.8 0.9 0.95"
MAF=0.01
THIN=100
```

and run it:

```bash
bash scripts/filter_vcf.sh
```

The script makes one VCF per missing-data threshold, so you can compare
them, and prints a table:

```
VCF                                  SNPs
my_data.vcf                          1270433
my_data-miss5-maf1-th100.recode.vcf  ...
my_data-miss8-maf1-th100.recode.vcf  ...
my_data-miss9-maf1-th100.recode.vcf  13936
my_data-miss95-maf1-th100.recode.vcf ...
```

(Example: 135 individuals, from 1.27 million raw variants to about 14,000
SNPs with `--max-missing 0.9`, MAF 1% and one SNP per 100 bp.)

**Which threshold to use?** There is a trade-off: stricter filters give
fewer SNPs but less missing data. A good practice is to run the main
analysis (e.g. PCA) with two or three thresholds. If the result is the same,
the filtering does not drive the result. `--max-missing 0.9` is a common
choice for RAD-seq.

Check the missing data per individual in the filtered file. Individuals with
much more missing data than the others can appear as separate groups in PCA
and STRUCTURE:

```bash
vcftools --vcf filtered/my_data-miss9-maf1-th100.recode.vcf --missing-indv --out my_data-miss9
sort -k5,5gr my_data-miss9.imiss | head
```

## 4. Optional: outlier SNPs with pcadapt

> **We no longer remove pcadapt outliers before the structure and
> phylogenetic analyses.** pcadapt flags SNPs that are strongly associated
> with the main PCs, that is, with population structure itself. Many of
> these are not under selection: they are simply the SNPs that differ most
> between groups. Removing them removes part of the signal you want to
> study, and can make groups look more similar than they are. Use the
> filtered VCF of section 3 for the other tutorials.
>
> pcadapt is still useful to **find candidate SNPs** for selection, or to
> check whether the results change without them. Report it as a
> sensitivity analysis, not as a filter.

The steps below show how to run it if you need the outlier list.

### 4.1 Choose K (number of PCs)

In `pcadapt_outliers.R`, set the input files and leave `K_FINAL <- NA`:

```r
VCF      <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
POPMAP   <- "popmap.txt"
K_MAX    <- 20
K_FINAL  <- NA
```

```bash
Rscript scripts/pcadapt_outliers.R
```

Open `pcadapt/my_data_pcadapt_screeplot.svg`. It shows the proportion of
variance explained by each PC. Choose **K = the number of PCs before the
curve becomes flat** (Cattell's rule): the PCs after that point describe
noise, not population structure.

Also look at the `_scores_PC*.svg` plots: if PC5 and PC6 still separate
groups of individuals, those PCs describe real structure and should be
included in K. If they show only a cloud of points, they do not. In the
example data the scree plot flattened around PC9, but the score plots
separated groups of individuals only up to PC5, so K = 5 was used.

### 4.2 Find the outliers

Set `K_FINAL` and run again:

```r
K_FINAL  <- 5
```

```bash
Rscript scripts/pcadapt_outliers.R
```

```
K = 5, bonferroni correction, alpha = 0.05: 757 outlier SNPs out of 13936 (5.4%)
Outlier positions: pcadapt/my_data_pcadapt_outliers.txt
```

The script writes the chromosome and position of each outlier, so there is
no need to convert the pcadapt SNP numbers to positions by hand.

pcadapt prints `Converter vcf to pcadapt is deprecated. Please use PLINK
for conversion to bed`. With a VCF of biallelic SNPs the conversion works:
the warning can be ignored. The script stops if pcadapt and the VCF do not
have the same number of SNPs.

- **Bonferroni** is strict and removes only strong outliers. `CORRECTION <-
  "BH"` (false discovery rate) removes more SNPs.
- Check the histogram of p-values (`_pvalues.svg`): it should be mostly
  flat, with a peak near 0 (the outliers). A U-shaped or very irregular
  histogram suggests that K is wrong.

### 4.3 Optional: a VCF without the outliers

```bash
vcftools --vcf filtered/my_data-miss9-maf1-th100.recode.vcf \
    --exclude-positions pcadapt/my_data_pcadapt_outliers.txt \
    --recode --recode-INFO-all --out filtered/my_data-miss9-maf1-th100.neutral
```

```
After filtering, kept 135 out of 135 Individuals
After filtering, kept 13179 out of a possible 13936 Sites
```

> **Do the outlier test separately for each VCF.** With or without the
> outgroup, the PCs are different, so the outliers are different too (in the
> example data: 757 outliers without the outgroup and 624 with it).

## 5. Check the result

- Keep a table of the filters and the number of samples and SNPs of each VCF
  for the methods section of your paper.

## References

- Danecek P, Auton A, Abecasis G, et al. (2011) The variant call format and
  VCFtools. *Bioinformatics* 27: 2156–2158.
- Danecek P, Bonfield JK, Liddle J, et al. (2021) Twelve years of SAMtools
  and BCFtools. *GigaScience* 10: giab008.
- Privé F, Luu K, Vilhjálmsson BJ, Blum MGB (2020) Performing highly
  efficient genome scans for local adaptation with R package pcadapt version
  4. *Molecular Biology and Evolution* 37: 2153–2154.
- Luu K, Bazin E, Blum MGB (2017) pcadapt: an R package to perform genome
  scans for selection based on principal component analysis. *Molecular
  Ecology Resources* 17: 67–77.
