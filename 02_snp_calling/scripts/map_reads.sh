#!/bin/bash
# Map the clean reads of each sample to the reference genome with BWA,
# keep only the mapped reads, sort and index the BAM files.
#
# Run from the analysis folder (02_snp_calling/):
#   conda activate snps
#   bash scripts/map_reads.sh

# ---- settings: edit these lines ----
REFERENCE=data/reference.fasta              # reference genome (indexed with bwa index, see the tutorial)
READS_DIR=../01_raw_reads/results/clean_reads   # clean reads (<sample>.fq.gz)
SAMPLES=data/samples.txt                    # list of sample names, one per line
THREADS=8

# ---- run ----
mkdir -p results/bam results/mapping_stats

for sample in $(cat $SAMPLES); do
    echo "Mapping $sample"
    # map | keep mapped reads (-F 4) | sort by position
    bwa mem -t $THREADS $REFERENCE $READS_DIR/${sample}.fq.gz \
        | samtools view -F 4 -b - \
        | samtools sort -@ $THREADS -o results/bam/${sample}.bam -
    samtools index results/bam/${sample}.bam
    samtools flagstat results/bam/${sample}.bam > results/mapping_stats/${sample}.flagstat.txt
done

echo "Done. BAM files are in results/bam/"
