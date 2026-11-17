#!/bin/bash
# Filter a VCF: keep the selected samples, then keep biallelic SNPs with
# little missing data, a minimum allele frequency and one SNP per locus.
#
# Run from the analysis folder (03_filtering/):
#   conda activate popgen
#   bash scripts/filter_vcf.sh

# ---- settings: edit these lines ----
INPUT=data/my_data.vcf          # VCF from 02_snp_calling
SAMPLES=data/samples_keep.txt   # samples to keep, one per line
PREFIX=my_data                  # name of the output files
MISSING=0.9                     # a site must be genotyped in at least 90% of the samples
MAF=0.01                        # minimum minor allele frequency
THIN=100                        # keep at most one SNP every 100 bp

# ---- 1. keep the selected samples and remove sites that became monomorphic ----
bcftools view -S $SAMPLES $INPUT -o results/${PREFIX}.subset.vcf
bcftools view -e 'AC==0 || AC==AN' results/${PREFIX}.subset.vcf -o results/${PREFIX}.vcf
rm results/${PREFIX}.subset.vcf

# ---- 2. filter ----
vcftools --vcf results/${PREFIX}.vcf \
    --min-alleles 2 --max-alleles 2 --remove-indels \
    --max-missing $MISSING --maf $MAF --thin $THIN \
    --recode --recode-INFO-all \
    --out results/${PREFIX}-miss${MISSING}-maf${MAF}-th${THIN}

# ---- 3. count ----
echo "Samples:              $(bcftools query -l results/${PREFIX}.vcf | wc -l)"
echo "Variants before:      $(grep -vc '^#' results/${PREFIX}.vcf)"
echo "SNPs after filtering: $(grep -vc '^#' results/${PREFIX}-miss${MISSING}-maf${MAF}-th${THIN}.recode.vcf)"
