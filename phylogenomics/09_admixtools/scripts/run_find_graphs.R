# Run find_graphs (ADMIXTOOLS 2) many times for each number of admixture
# events, in parallel, and save the best graph of each run.
#
# Usage (in the results/ folder, inside byobu, with the conda env active):
#   Rscript ../scripts/run_find_graphs.R
#
# Output: one file per run, winner.<replicate>m<numadmix>.rds
#   (e.g. winner.12m3.rds = replicate 12 with 3 admixture events)

library(admixtools)
library(dplyr)
library(doSNOW)
library(foreach)

# ---- Settings: edit here --------------------------------------------------
my_f2_dir    <- "f2_my_data"   # folder made by extract_f2()
outpop       <- "Outgroup"     # outgroup population (as in the .fam file)
numadmix_all <- 0:6            # numbers of admixture events to test
replicates   <- 1:100          # number of independent runs per numadmix
stop_gen2    <- 1000           # maximum number of generations per run
numgraphs    <- 100            # graphs per generation
cores        <- 12             # number of runs at the same time
# ---------------------------------------------------------------------------

cl <- makeCluster(cores, outfile = "")
registerDoSNOW(cl)

run_find_graphs <- function(my_f2_dir, outpop, numadmix, replicate_num) {
  out_file <- paste0("winner.", replicate_num, "m", numadmix, ".rds")
  if (file.exists(out_file)) return(NULL)    # already done: skip (useful to restart)

  print(paste("Starting numadmix =", numadmix, "replicate =", replicate_num))
  high <- admixtools::find_graphs(
    my_f2_dir,
    outpop = outpop,
    numadmix = numadmix,
    stop_gen2 = stop_gen2,
    numgraphs = numgraphs,
    plusminus_generations = 20,
    verbose = FALSE
  )
  # keep the graph with the lowest score of this run
  winner <- high %>% dplyr::slice_min(score, with_ties = FALSE)
  saveRDS(winner, file = out_file)
  print(paste("Completed numadmix =", numadmix, "replicate =", replicate_num))
  NULL
}

foreach(numadmix = numadmix_all) %:%
  foreach(replicate_num = replicates,
          .export = c("stop_gen2", "numgraphs")) %dopar% {
    run_find_graphs(my_f2_dir, outpop, numadmix, replicate_num)
  }

stopCluster(cl)
cat("Done\n")
