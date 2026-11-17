#!/bin/bash
# Merge the BAM files of all samples (adding the sample name to the reads)
# and call the variants with FreeBayes.
#
# Run from the analysis folder (02_snp_calling/):
#   conda activate snps
#   bash scripts/call_snps.sh

# ---- settings: edit these lines ----
REFERENCE=data/reference.fasta
SAMPLES=data/samples.txt            # list of sample names, one per line (same as map_reads.sh)
OUTPUT=results/my_data.vcf          # VCF with all variants
MIN_MAP_QUALITY=30                  # -m: minimum mapping quality of a read
MIN_BASE_QUALITY=30                 # -q: minimum quality of a base
MIN_COVERAGE=10                     # --min-coverage: minimum depth at a site (all samples together)
LIMIT_COVERAGE=10000                # --limit-coverage: skip sites with higher depth (repeats)

# ---- 1. merge the BAM files ----
# bamaddrg needs "-b <file> -s <sample name>" for each sample
BAM_LIST=""
for sample in $(cat $SAMPLES); do
    BAM_LIST="$BAM_LIST -b results/bam/${sample}.bam -s ${sample}"
done

echo "Merging BAM files"
bamaddrg $BAM_LIST > results/all_samples.bam
samtools index results/all_samples.bam
samtools faidx $REFERENCE

# ---- 2. call the variants ----
echo "Calling variants with FreeBayes (this can take hours)"
freebayes -b results/all_samples.bam -f $REFERENCE -v $OUTPUT \
    -m $MIN_MAP_QUALITY -q $MIN_BASE_QUALITY \
    --min-coverage $MIN_COVERAGE --limit-coverage $LIMIT_COVERAGE

# ---- 3. summary ----
echo "Number of variants: $(grep -vc '^#' $OUTPUT)"
echo "Number of samples:  $(bcftools query -l $OUTPUT | wc -l)"
bcftools stats $OUTPUT > results/my_data.vcf.stats.txt
echo "Done: $OUTPUT"
