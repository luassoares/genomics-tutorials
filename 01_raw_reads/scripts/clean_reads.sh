#!/bin/bash
# Remove barcodes and adapters, trim low-quality bases and discard short
# reads with fastq-mcf, one sample at a time.
#
# Run from the analysis folder (01_raw_reads/):
#   conda activate reads
#   bash scripts/clean_reads.sh

# ---- settings: edit these lines ----
RAW_DIR=../00_raw_data                  # folder with the raw reads (<sample>.FASTQ.gz)
SAMPLES=data/samples.txt                # list of sample names, one per line (file name without .FASTQ.gz)
ADAPTERS=data/barcodes_adapter.fasta    # barcodes + Illumina adapter (see the tutorial)
QUALITY=30                              # trim bases with quality below this value
MIN_LENGTH=50                           # discard reads shorter than this after trimming

# ---- run ----
mkdir -p results/clean_reads results/logs

for sample in $(cat $SAMPLES); do
    echo "Cleaning $sample"
    fastq-mcf -q $QUALITY -l $MIN_LENGTH \
        -o results/clean_reads/${sample}.fq.gz \
        $ADAPTERS $RAW_DIR/${sample}.FASTQ.gz > results/logs/${sample}.fastq-mcf.log
done

echo "Done. Clean reads are in results/clean_reads/"
