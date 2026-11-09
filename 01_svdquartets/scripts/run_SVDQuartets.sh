#!/bin/bash
# SVDQuartets (PAUP*): best tree + bootstrap consensus, rooted on an outgroup.
#
# Usage: edit the settings below, then
#        bash run_SVDQuartets.sh        (run inside byobu/tmux/screen for long runs)
#
# Output (in ${OUTDIR}/):
#   ${PREFIX}.nex              PAUP* batch file written by this script
#   ${PREFIX}.log              PAUP* log
#   ${PREFIX}.out              screen output, with start and end times
#   ${PREFIX}_besttree.tre     best tree (fully resolved)
#   ${PREFIX}_bootstrap.tre    bootstrap trees
#   ${PREFIX}_consensus.tre    majority-rule consensus (PAUP* does not write the
#                              support values: use add_support_and_pops.py)

# ---- settings ----
INPUT=my_data.Ind.nexus          # NEXUS alignment (taxon names must not be only numbers)
PREFIX=SVDQ_my_data              # prefix of the output files
NREPS=1000                       # bootstrap replicates
NTHREADS=8                       # CPU threads (check with nproc)
OUTDIR=results
CONDA_ENV=paup                   # conda environment with libgfortran4 (empty = do not activate)

# Outgroup: either list the taxon names...
OUTGROUP=""                      # e.g. "Ind_101 Ind_102 Ind_103"
# ...or take them from a popmap (<population> <tab> <sample ID>):
POPMAP=""                        # e.g. ../popmap.txt
OUTGROUP_POP=""                  # population name of the outgroup in the popmap
ID_PREFIX="Ind_"                 # prefix added to the sample IDs in the NEXUS file

# Species tree instead of a lineage tree: name of a taxpartition defined in a
# SETS block of the input file (empty = individuals as tips)
TAXPARTITION=""

# ---- activate conda env ----
# (conda activate does not work in a script without sourcing conda.sh first)
if [ -n "${CONDA_ENV}" ]; then
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate ${CONDA_ENV}
fi

# ---- checks ----
[ -f "${INPUT}" ] || { echo "Input file not found: ${INPUT}"; exit 1; }
command -v paup > /dev/null || { echo "paup not found in PATH (see 00_setup)"; exit 1; }

if [ -z "${OUTGROUP}" ] && [ -n "${POPMAP}" ] && [ -n "${OUTGROUP_POP}" ]; then
    # tr removes Windows line endings (CRLF) from the popmap
    OUTGROUP=$(tr -d '\r' < ${POPMAP} | awk -v p="${OUTGROUP_POP}" -v pre="${ID_PREFIX}" '$1 == p {printf "%s%s ", pre, $2}')
fi
[ -n "${OUTGROUP}" ] || { echo "No outgroup defined"; exit 1; }

PARTITION_OPT=""
[ -n "${TAXPARTITION}" ] && PARTITION_OPT="taxpartition=${TAXPARTITION}"

# absolute paths, so PAUP* finds the files wherever it runs from
mkdir -p ${OUTDIR}
INPUT=$(realpath ${INPUT})
OUT=$(realpath ${OUTDIR})/${PREFIX}

# ---- write PAUP* batch file ----
cat > ${OUT}.nex << EOF
#NEXUS
begin paup;
  set autoclose=yes warntree=no warnreset=no increase=auto;
  log start file=${OUT}.log replace;
  execute ${INPUT};

  outgroup ${OUTGROUP};

  [1. best tree (fully resolved, no polytomies)]
  svdquartets evalQuartets=all ambigs=distribute nthreads=${NTHREADS} ${PARTITION_OPT};
  roottrees;
  savetrees file=${OUT}_besttree.tre format=nexus root=yes brLens=no replace;

  [2. bootstrap - majority-rule consensus]
  svdquartets evalQuartets=all ambigs=distribute bootstrap=standard nreps=${NREPS} nthreads=${NTHREADS} ${PARTITION_OPT}
              treeFile=${OUT}_bootstrap.tre replace;
  roottrees;
  savetrees file=${OUT}_consensus.tre format=nexus root=yes brLens=no supportValues=nodeLabels replace;

  log stop;
  quit;
end;
EOF

# ---- run (show on screen and save to ${OUT}.out) ----
exec > >(tee ${OUT}.out) 2>&1
echo "Start: $(date)"
echo "Outgroup: ${OUTGROUP}"
paup -n ${OUT}.nex
echo "End: $(date)"
