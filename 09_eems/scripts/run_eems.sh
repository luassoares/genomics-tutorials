#!/bin/bash
# Run EEMS (runeems_snps) with several independent chains: writes one params
# file per chain and runs the chains in parallel.
#
# Usage: edit the settings below, then
#        bash run_eems.sh          (inside byobu/tmux/screen: it takes hours)
#
# Input: ${DATAPATH}.diffs, ${DATAPATH}.coord (from vcf2diffs.R), ${DATAPATH}.outer
# Output (in ${OUTDIR}/):
#   params-chain<N>.ini           EEMS parameter files
#   ${NAME}-nDemes${NDEMES}-chain<N>/   EEMS results of each chain (input of plot_eems.R)
#   ${NAME}-chain<N>.log          screen output of each chain

# ---- settings ----
DATAPATH=eems/data/my_data           # same as in vcf2diffs.R (without extension)
NAME=my_data                         # prefix of the output folders
OUTDIR=eems/runs
NINDIV=131                           # printed by vcf2diffs.R
NSITES=13179                         # printed by vcf2diffs.R
NDEMES=300                           # density of the grid (approximate number of demes)
NCHAINS=3                            # independent chains, run in parallel
SEEDS="123 456 789"                  # one seed per chain
NUM_MCMC=4000000                     # MCMC iterations
NUM_BURN=1000000                     # burn-in iterations
NUM_THIN=9999                        # thinning interval
EEMS=runeems_snps                    # path of the runeems_snps binary, if not in PATH

# ---- checks ----
for ext in diffs coord outer; do
    [ -f "${DATAPATH}.${ext}" ] || { echo "Missing ${DATAPATH}.${ext}"; exit 1; }
done
command -v ${EEMS} > /dev/null || [ -x "${EEMS}" ] || { echo "${EEMS} not found (see section 2 of the tutorial)"; exit 1; }
n=$(wc -l < ${DATAPATH}.coord)
[ "${n}" -eq "${NINDIV}" ] || { echo "NINDIV=${NINDIV} but ${DATAPATH}.coord has ${n} lines"; exit 1; }
read -ra SEED_ARR <<< "${SEEDS}"
[ "${#SEED_ARR[@]}" -ge "${NCHAINS}" ] || { echo "Give at least ${NCHAINS} seeds"; exit 1; }

mkdir -p ${OUTDIR}
DATAPATH=$(realpath -m ${DATAPATH})     # absolute paths: EEMS runs from any folder
OUTDIR=$(realpath ${OUTDIR})

# ---- write params and run ----
for i in $(seq 1 ${NCHAINS}); do
    PARAMS=${OUTDIR}/params-chain${i}.ini
    cat > ${PARAMS} << EOF
datapath = ${DATAPATH}
mcmcpath = ${OUTDIR}/${NAME}-nDemes${NDEMES}-chain${i}
nIndiv = ${NINDIV}
nSites = ${NSITES}
nDemes = ${NDEMES}
diploid = true
numMCMCIter = ${NUM_MCMC}
numBurnIter = ${NUM_BURN}
numThinIter = ${NUM_THIN}
EOF
    echo "Chain ${i}: seed ${SEED_ARR[$((i-1))]}, start $(date)"
    ${EEMS} --params ${PARAMS} --seed ${SEED_ARR[$((i-1))]} > ${OUTDIR}/${NAME}-chain${i}.log 2>&1 &
done

wait
echo "All chains finished: $(date)"
for i in $(seq 1 ${NCHAINS}); do
    tail -n 3 ${OUTDIR}/${NAME}-chain${i}.log
done
