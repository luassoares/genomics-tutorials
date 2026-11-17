# Tutorials: population genomics and phylogenomics with SNP data

Step-by-step tutorials, with simple scripts, for analysing RAD-seq or
DArTseq data of closely related populations and species: from the raw
reads to the population structure, the relationships between groups and
the tests of gene flow. They are written for students, on a Linux server
without administrator (root) access, with conda, and RStudio for the
analyses in R.

## Part 1 — From raw reads to population structure

Follow the tutorials in order: the output of each one is the input of the
next.

| # | Tutorial | Question |
|---|---|---|
| 00 | [Linux basics](00_linux_basics/) | How to use the terminal: main commands, byobu, the PATH |
| 00 | [Setup](00_setup/) | How to install the programs (conda) and organize the project folders |
| 01 | [Raw reads](01_raw_reads/) | Are the reads good? Remove barcodes, adapters and low-quality bases |
| 02 | [Mapping and SNP calling](02_snp_calling/) | From clean reads to a VCF (reference genome, or Stacks de novo) |
| 03 | [Filtering](03_filtering/) | Which samples and SNPs to keep? |
| 04 | [Genetic diversity](04_genetic_diversity/) | How diverse is each population, and how different are they (π, H<sub>o</sub>, F<sub>IS</sub>, F<sub>ST</sub>)? |
| 05 | [PCA and DAPC](05_pca_dapc/) | How are the individuals grouped, without a model? (RStudio) |
| 06 | [STRUCTURE and fastSTRUCTURE](06_structure/) | How many genetic clusters, and what is the ancestry of each individual? (plots with pophelper) |
| 07 | [sNMF and TESS3](07_snmf_tess3/) | Ancestry coefficients, with and without geography, on a map |
| 08 | [conStruct](08_construct/) | Discrete groups, isolation by distance, or both? |
| 09 | [EEMS](09_eems/) | Where are the barriers and corridors to gene flow? |
| 10 | [FEEMS](10_feems/) | The same, faster, with cross-validation |

## Part 2 — Phylogenomics and introgression

In the folder [phylogenomics](phylogenomics/):

| # | Tutorial | Question |
|---|---|---|
| 01 | [SVDQuartets](phylogenomics/01_svdquartets/) | How are the individuals or populations related? |
| 02 | [Ancestry painting](phylogenomics/02_ancestry_painting/) | Are putative hybrids heterozygous at sites fixed between the parents? |
| 03 | [Dsuite](phylogenomics/03_dsuite/) | Is there excess allele sharing (gene flow) between groups? |

More tutorials will be added.

## How to use

1. New to the terminal? Start with [00_linux_basics](00_linux_basics/).
   Then read [00_setup](00_setup/): install conda and create the project
   folders.
2. Each tutorial starts with an **Install** section (the programs needed for
   that analysis only) and **Prepare the folder** (which input files go in
   `data/`).
3. Copy the scripts of the tutorial to the `scripts/` folder of your
   analysis, change the settings at the top (file names, number of
   processors), and run them as shown in the tutorial. The R scripts are
   meant to be run in RStudio, line by line.

```bash
cd ~
git clone https://github.com/luassoares/genomics-tutorials.git
```

The examples come from analyses of *Petunia* species (DArTseq data, ~50–150
individuals, ~10,000 SNPs after filtering). Each tutorial ends with links to
more detailed tutorials and the references to cite.

## Credits

The phylogenomics tutorials adapt and extend the tutorials of **Milan
Malinsky** and **Michael Matschiner** (links in
[phylogenomics](phylogenomics/)). Please cite the original methods, listed
at the end of each tutorial.

## Author

Luana Soares
