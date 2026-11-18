# Part 2 — Phylogenomics and introgression

Tutorials to infer the relationships between individuals, populations and
species, estimate divergence times, and test for hybridization and gene
flow. They start from the filtered VCF of [03_filtering](../03_filtering/)
and use the results of Part 1 (which groups exist) to define the tips and
the hypotheses.

| # | Tutorial | Question |
|---|---|---|
| 01 | [SplitsTree](01_splitstree/) | Is the signal tree-like, or is there conflict (reticulation)? A first look with a network |
| 02 | [RAxML-ng](02_raxml/) | Maximum-likelihood tree of the individuals (concatenated SNPs) |
| 03 | [SVDQuartets](03_svdquartets/) | Lineage and species tree under the multispecies coalescent |
| 04 | [SNAPP](04_snapp/) | Species tree with divergence times |
| 05 | [Ancestry painting](05_ancestry_painting/) | Are putative hybrids heterozygous at sites fixed between the parents? |
| 06 | [Dsuite](06_dsuite/) | Is there excess allele sharing (gene flow) between groups? D statistics, f4-ratio, f-branch |

More tutorials will be added: HyDe, Dfoil and ADMIXTOOLS (introgression),
SNaQ and PhyloNet (phylogenetic networks), IQ-TREE with concordance factors,
and species delimitation (BPP, DELINEATE, HHSD, delimitR).

## Before starting

- **Install and folders**: as in Part 1 ([00_setup](../00_setup/)). Each
  tutorial has its own **Install** section. Keep the analyses in a
  `phylogenomics/` folder of your project, one folder per analysis with
  `data/`, `scripts/` and `results/`.
- **An outgroup** is needed to root the trees and for the D statistics. Make
  a VCF **with the outgroup** samples in [03_filtering](../03_filtering/#1-which-vcf-for-which-analysis):
  add them to `samples_keep.txt` and use another prefix, e.g.
  `PREFIX=my_data-out` → `my_data-out-miss0.9-maf0.01-th100.recode.vcf`.
- **Thinned or not**: the tree methods use the thinned VCF (one SNP per
  locus); ancestry painting and Dsuite use a VCF that is not thinned (see
  each tutorial).

## More information

These tutorials adapt and extend the tutorials of **Milan Malinsky** and
**Michael Matschiner**:

- Species tree inference with SNP data:
  <https://github.com/millanek/tutorials/tree/master/species_tree_inference_with_snp_data>
- Analysis of introgression with SNP data:
  <https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data>
- Divergence time estimation with SNP data:
  <https://github.com/mmatschiner/tutorials/blob/master/divergence_time_estimation_with_snp_data/README.md>
- Scripts from those tutorials (`get_fixed_site_gts.rb`,
  `plot_fixed_site_gts.rb`, `snapp_prep.rb`), `vcf2phylip.py` and Dsuite's
  `dtools.py` are downloaded from their original repositories and are not
  included here.
- Review of methods to detect introgression: Hibbins MS, Hahn MW (2022)
  Phylogenomic approaches to detecting and characterizing introgression.
  *Genetics* 220: iyab173.
