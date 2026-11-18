# Run SNaQ for h = 0, 1, 2, ... hybridizations.
# Each network starts from the best network with one hybridization less.
#
# Usage (in the results/ folder, inside byobu, with the conda env active):
#   julia ../scripts/run_snaq.jl
#
# Output, in snaq/: net<h>.out (best network and its score), net<h>.log,
#   net<h>.networks (other directions of the hybrid edges), net<h>.err
# If a run stops, you can start again from the last net<h>.out
# (set hmin and start_net below).

using Distributed

# ---- Settings: edit here --------------------------------------------------
cf_file    = "my_data_CF.csv"       # table from SNPs2CF
start_tree = "../data/start_tree.nwk" # population tree (tips = species in the Imap)
hmin       = 0                      # first h
hmax_all   = 6                      # last h
nruns      = 10                     # independent runs for each h (default 10)
ncores     = 10                     # processors (one run per processor)
start_net  = ""                     # to restart: e.g. "snaq/net3.out" and hmin = 4
# ---------------------------------------------------------------------------

addprocs(ncores)
@everywhere using PhyloNetworks, SNaQ

mkpath("snaq")
cf = readtableCF(cf_file)

if start_net == ""
    net = readnewick(start_tree)
else
    net = readsnaqnetwork(start_net)
end

for h in hmin:hmax_all
    println("\n========== h = $h ==========")
    global net = snaq!(net, cf, hmax=h, filename="snaq/net$h",
                       seed=1000 + h, runs=nruns)
end
println("Done")
