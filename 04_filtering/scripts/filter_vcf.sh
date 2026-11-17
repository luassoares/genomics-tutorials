#!/bin/bash
# Subset a VCF to a list of samples and filter it with several missing-data
# thresholds at once; writes a table with the number of SNPs kept per file.
#
# Usage: edit the settings below, then
#        bash filter_vcf.sh
#
# Output (in ${OUTDIR}/):
#   ${PREFIX}.vcf                              VCF with the selected samples, monomorphic sites removed
#   ${PREFIX}-miss<M>-maf<F>-th<T>.recode.vcf  one filtered VCF per missing-data threshold
#   ${PREFIX}_snp_counts.txt                   number of SNPs in each file
#   ${PREFIX}_filtering.out                    screen output

# ---- settings ----
INPUT=all_samples.vcf            # VCF from the SNP caller (freebayes, Stacks...)
PREFIX=my_data                   # prefix of the output files
OUTDIR=filtered

SAMPLES=""                       # file with one sample ID per line to KEEP (empty = all samples)
EXCLUDE=""                       # comma-separated sample IDs to REMOVE, e.g. duplicates: "1023,1104"

MISSING="0.5 0.8 0.9 0.95"       # --max-missing values: 0.9 = a site must be genotyped in 90% of the samples
MAF=0.01                         # minor allele frequency (0.01 = 1%)
THIN=100                         # keep at most one SNP every THIN bp (0 = no thinning)

# ---- checks ----
[ -f "${INPUT}" ] || { echo "Input file not found: ${INPUT}"; exit 1; }
for tool in bcftools vcftools; do
    command -v ${tool} > /dev/null || { echo "${tool} not found in PATH (see 00_setup)"; exit 1; }
done
[ -z "${SAMPLES}" ] || [ -f "${SAMPLES}" ] || { echo "Sample list not found: ${SAMPLES}"; exit 1; }

mkdir -p ${OUTDIR}
exec > >(tee ${OUTDIR}/${PREFIX}_filtering.out) 2>&1
echo "Start: $(date)"

# ---- 1. subset samples and remove sites that became monomorphic ----
# After removing samples, some sites have only one allele left (AC==0: only
# reference; AC==AN: only alternative). They carry no information.
SUBSET_OPT=""
if [ -n "${SAMPLES}" ]; then
    tr -d '\r' < ${SAMPLES} > ${OUTDIR}/.samples_unix.txt      # remove Windows line endings
    SUBSET_OPT="-S ${OUTDIR}/.samples_unix.txt"
fi
[ -n "${EXCLUDE}" ] && SUBSET_OPT="${SUBSET_OPT} -s ^${EXCLUDE}"

bcftools view ${SUBSET_OPT} ${INPUT} -Ou \
    | bcftools view -e 'AC==0 || AC==AN' -o ${OUTDIR}/${PREFIX}.vcf
rm -f ${OUTDIR}/.samples_unix.txt

NSAMPLES=$(bcftools query -l ${OUTDIR}/${PREFIX}.vcf | wc -l)
echo "Samples kept: ${NSAMPLES}"

# ---- 2. filter: biallelic SNPs, missing data, MAF, thinning ----
# file-name labels: 0.9 -> miss9, 0.95 -> miss95, 0.01 -> maf1
MAF_LABEL=$(echo ${MAF} | sed 's/^0\.0*//')
THIN_OPT=""
THIN_LABEL=""
if [ "${THIN}" -gt 0 ]; then
    THIN_OPT="--thin ${THIN}"
    THIN_LABEL="-th${THIN}"
fi

echo -e "VCF\tSNPs" > ${OUTDIR}/${PREFIX}_snp_counts.txt
echo -e "${PREFIX}.vcf\t$(grep -vc '^#' ${OUTDIR}/${PREFIX}.vcf)" >> ${OUTDIR}/${PREFIX}_snp_counts.txt

for miss in ${MISSING}; do
    OUT=${PREFIX}-miss${miss#0.}-maf${MAF_LABEL}${THIN_LABEL}
    vcftools --vcf ${OUTDIR}/${PREFIX}.vcf \
        --min-alleles 2 --max-alleles 2 --remove-indels \
        --max-missing ${miss} --maf ${MAF} ${THIN_OPT} \
        --recode --recode-INFO-all --out ${OUTDIR}/${OUT}
    echo -e "${OUT}.recode.vcf\t$(grep -vc '^#' ${OUTDIR}/${OUT}.recode.vcf)" >> ${OUTDIR}/${PREFIX}_snp_counts.txt
done

echo
column -t ${OUTDIR}/${PREFIX}_snp_counts.txt 2>/dev/null || cat ${OUTDIR}/${PREFIX}_snp_counts.txt
echo "End: $(date)"
