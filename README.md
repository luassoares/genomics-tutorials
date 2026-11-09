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
  the sample IDs exactly as in the VCF:

  ```
  PopA	1023
  PopA	1024
  PopB	1101
  Outgroup1	1201
  ```

## Suggested project structure

```
my_project/
├── 01_vcfs/
├── 01_popmap/
├── 02_SVDQuartets/
├── 03_AncestryPainting/
└── 04_Dsuite/
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
