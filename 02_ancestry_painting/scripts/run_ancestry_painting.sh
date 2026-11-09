#!/bin/bash
# Ancestry painting (scripts by M. Matschiner, from the introgression tutorial of
# M. Malinsky & M. Matschiner):
# https://github.com/millanek/tutorials/tree/master/analysis_of_introgression_with_snp_data
#
# Runs get_fixed_site_gts.rb and plot_fixed_site_gts.rb for one or more hybrid
# scenarios. Sample lists are taken from a popmap (<population> <tab> <sample ID>)
# with regular expressions on the population names.
#
# Usage: edit the settings below, then
#        bash run_ancestry_painting.sh
#        Output is shown on screen and saved to ${OUTDIR}/ancestry_painting.out

# ---- settings ----
VCF=my_data.biallelic.vcf         # uncompressed VCF (gunzip -c file.vcf.gz > file.vcf)
POPMAP=popmap.txt
OUTDIR=results
PARENT_COMPLETENESS=0.8   # min. proportion of parental alleles with data (1.0 = no missing data)
PARENT_FIXEDNESS=1.0      # min. frequency of the fixed allele in each parent (1.0 = completely fixed)
PLOT_COMPLETENESS=0.8     # min. completeness at a site, parents and hybrids, for the plot
MIN_DIST=100              # min. distance (bp) between plotted SNPs

P1_REGEX='^ParentA'       # populations of parent 1 (regular expression on column 1 of the popmap)
P2_REGEX='^ParentB'       # populations of parent 2
# hybrid scenarios, name:population-regex, separated by spaces
SCENARIOS="ALL_HYB:^Hyb HYB1:^Hyb1$ HYB2:^Hyb2$"

SRC=https://raw.githubusercontent.com/millanek/tutorials/master/analysis_of_introgression_with_snp_data/src

# ---- checks ----
command -v ruby > /dev/null || { echo "ruby not found (install with: conda install -c conda-forge ruby)"; exit 1; }
[ -f "${VCF}" ] || { echo "VCF not found: ${VCF}"; exit 1; }
[ -f "${POPMAP}" ] || { echo "Popmap not found: ${POPMAP}"; exit 1; }
for f in get_fixed_site_gts.rb plot_fixed_site_gts.rb; do
    [ -f $f ] || wget -q ${SRC}/$f
done
mkdir -p ${OUTDIR}

# ---- show on screen and save to log ----
exec > >(tee ${OUTDIR}/ancestry_painting.out) 2>&1
echo "Start: $(date)"
echo "VCF: ${VCF}"

# comma-separated sample IDs of the populations matching a regex (column 1 of the popmap)
ids() { tr -d '\r' < ${POPMAP} | awk -v re="$1" '$1 ~ re {print $2}' | paste -sd, -; }

P1=$(ids "${P1_REGEX}")
P2=$(ids "${P2_REGEX}")
[ -n "${P1}" ] && [ -n "${P2}" ] || { echo "No samples match ${P1_REGEX} or ${P2_REGEX} in ${POPMAP}"; exit 1; }
echo "Parent 1 (${P1_REGEX}): $(echo ${P1} | tr ',' '\n' | wc -l) individuals"
echo "Parent 2 (${P2_REGEX}): $(echo ${P2} | tr ',' '\n' | wc -l) individuals"

# ---- run each scenario ----
for sc in ${SCENARIOS}; do
    name=${sc%%:*}
    regex=${sc#*:}
    HYB=$(ids "${regex}")
    echo
    echo "==== ${name}: hybrids ${regex} ($(echo ${HYB} | tr ',' '\n' | wc -l) individuals) ===="

    ruby get_fixed_site_gts.rb ${VCF} ${OUTDIR}/${name}.fixed.txt \
        ${P1} ${HYB} ${P2} ${PARENT_COMPLETENESS} ${PARENT_FIXEDNESS}

    echo "Fixed sites: $(( $(wc -l < ${OUTDIR}/${name}.fixed.txt) - 1 ))"

    ruby plot_fixed_site_gts.rb ${OUTDIR}/${name}.fixed.txt ${OUTDIR}/${name}.fixed.svg \
        ${PLOT_COMPLETENESS} ${MIN_DIST}
done

echo
echo "End: $(date)"
