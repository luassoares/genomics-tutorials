# Tutorials: phylogenomics and introgression with SNP data

Step-by-step tutorials, with ready-to-use scripts, for analysing RAD-seq or
other reduced-representation SNP data (VCF) from closely related populations
and putative hybrids. They are written for a Linux server without
administrator (root) access.

| # | Tutorial | Question |
|---|---|---|
| 00 | [Setup](00_setup/) | Install PAUP\*, Dsuite and the other tools |
| 01 | [SVDQuartets](01_svdquartets/) | How are the individuals or populations related? |
| 02 | [Ancestry painting](02_ancestry_painting/) | Are putative hybrids heterozygous at sites fixed between the parents? |
| 03 | [Dsuite](03_dsuite/) | Is there excess allele sharing (gene flow) between groups? D statistics, f4-ratio, f-branch |
| 04 | [Filtering](04_filtering/) | Which samples and SNPs to keep? |
| 05 | [Genetic diversity](05_genetic_diversity/) | How diverse is each population, and how differentiated are they (π, H<sub>o</sub>, F<sub>IS</sub>, F<sub>ST</sub>)? |
| 06 | [PCA and DAPC](06_pca_dapc/) | How are individuals grouped, without a model? |
| 07 | [STRUCTURE and fastSTRUCTURE](07_structure/) | How many genetic clusters, and what is the ancestry of each individual? |
| 08 | [conStruct](08_construct/) | Discrete clusters, isolation by distance, or both? |
| 09 | [EEMS](09_eems/) | Where are the barriers and corridors to gene flow? |
| 10 | [sNMF and TESS3](10_snmf_tess3/) | Ancestry coefficients, with and without geography, on a map |
| 11 | [FEEMS](11_feems/) | Barriers and corridors to gene flow, fast, with cross-validation |

### Suggested order

The numbers follow the order in which the tutorials were written. For a new
dataset, a typical order is:

1. **Prepare the data:** [00 Setup](00_setup/) → [04 Filtering](04_filtering/)
2. **Describe the populations:** [05 Diversity](05_genetic_diversity/) →
   [06 PCA/DAPC](06_pca_dapc/) → [07 STRUCTURE](07_structure/) →
   [10 sNMF/TESS3](10_snmf_tess3/) → [08 conStruct](08_construct/) →
   [09 EEMS](09_eems/) → [11 FEEMS](11_feems/)
3. **Relationships:** [01 SVDQuartets](01_svdquartets/)
4. **Hybridization and gene flow:** [02 Ancestry painting](02_ancestry_painting/) →
   [03 Dsuite](03_dsuite/)

More tutorials will be added.

## How to use

Each folder has a `README.md` (the tutorial) and a `scripts/` folder. The
scripts have their settings at the top: copy the script to your analysis
folder, edit the settings and run it.

```bash
git clone https://github.com/<user>/<repository>.git
```

### Input files used in all tutorials

- **VCF** with the SNPs of all individuals.
- **Popmap**: two columns, tab-separated, `<population> <sample ID>`, with
  the sample IDs exactly as in the VCF (Stacks uses the opposite order, see
  [05](05_genetic_diversity/)):

  ```
  PopA	1023
  PopA	1024
  PopB	1101
  Outgroup1	1201
  ```

- **Coordinates** (for conStruct, EEMS, TESS3 and FEEMS): `<sample ID> <tab> <longitude>
  <tab> <latitude>`, decimal degrees, longitude first:

  ```
  1023	-50.950	-25.200
  1101	-51.217	-25.333
  ```

## Suggested project structure

```
my_project/
├── 01_vcfs/              raw and filtered VCFs (04)
├── 01_popmap/            popmap, coordinates
├── 02_SVDQuartets/
├── 03_AncestryPainting/
├── 04_Dsuite/
├── 05_diversity/
├── 06_pca_dapc/
├── 07_structure/
├── 08_construct/
├── 09_eems/
├── 10_snmf_tess3/
└── 11_feems/
```

## Credits

These tutorials adapt and extend the tutorials of **Milan Malinsky** and
**Michael Matschiner**. Please cite the original methods (listed in each
tutorial) and see the original tutorials:

- Species tree inference with SNP data:
  <https://github.com/millanek/tutorials/tree/master/species_tree_inference_with_snp_data>
- Analysis of introgression with SNP data:
  <https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data>
- Divergence time estimation with SNP data:
  <https://github.com/millanek/tutorials/tree/master/divergence_time_estimation_with_snp_data>

Scripts from those tutorials (`vcf2phylip.py`, `get_fixed_site_gts.rb`,
`plot_fixed_site_gts.rb`) and from Dsuite (`dtools.py`) are downloaded from
their original repositories and are not included here.

## Author

Luana Soares
