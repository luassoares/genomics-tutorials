#!/usr/bin/env Rscript
# Summarise and plot STRUCTURE or fastSTRUCTURE results (from structure_threader)
# with pophelper: Evanno table (STRUCTURE), replicates aligned and merged per K,
# and bar plots with population labels.
#
# Usage: edit the settings below, then
#   Rscript plot_structure.R
#
# Output (in OUTDIR, file names start with PREFIX):
#   _summary.txt         number of runs, mean likelihood (STRUCTURE) per K
#   _evanno.txt          Evanno table (delta K) — STRUCTURE only, needs >= 3 values of K
#   evannoMethodStructure.pdf  Evanno plots — STRUCTURE only
#   _K<min>-<max>.pdf    bar plots of the merged replicates, one row per K

# ---- settings ----
RESULTS_DIR <- "structure_out"              # output folder of structure_threader (-o)
PROGRAM     <- "structure"                  # "structure" or "faststructure"
INDIVIDUALS <- "my_data_individuals.txt"    # from vcf2structure.py: <ID> <tab> <population>, same order as the .str
PREFIX      <- "my_data_structure"
OUTDIR      <- "structure_plots"
K_PLOT      <- 2:6                          # values of K to plot
COLORS      <- NULL                         # NULL = pophelper default, or e.g. c("#5c9ead", "#da8e95", ...)

# ---- packages ----
if (!requireNamespace("pophelper", quietly = TRUE)) {
  stop("Install pophelper: remotes::install_github('royfrancis/pophelper') (see 00_setup)")
}
suppressPackageStartupMessages(library(pophelper))
dir.create(OUTDIR, showWarnings = FALSE)
out <- function(x) file.path(OUTDIR, paste0(PREFIX, x))

# ---- read results ----
pattern <- if (PROGRAM == "structure") "_f$" else "\\.meanQ$"
files <- list.files(RESULTS_DIR, pattern = pattern, full.names = TRUE, recursive = TRUE)
if (length(files) == 0) stop("No result files matching '", pattern, "' in ", RESULTS_DIR)
message(length(files), " result files")

slist <- readQ(files = files)

inds <- read.table(INDIVIDUALS, sep = "\t", header = FALSE, colClasses = "character",
                   col.names = c("id", "pop"))
n_ind <- unique(sapply(slist, nrow))
if (length(n_ind) != 1 || n_ind != nrow(inds)) {
  stop("The runs have ", paste(n_ind, collapse = "/"), " individuals but ", INDIVIDUALS,
       " has ", nrow(inds), ". Use the file made together with the .str input.")
}
slist <- lapply(slist, "rownames<-", inds$id)

# ---- summary and Evanno ----
tr <- tabulateQ(slist)
sr <- summariseQ(tr)
write.table(sr, out("_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sr)

if (PROGRAM == "structure") {
  ev <- tryCatch(
    evannoMethodStructure(data = sr, exportplot = TRUE, exportpath = OUTDIR, imgtype = "pdf"),
    error = function(e) { message("Evanno method not possible: ", conditionMessage(e)); NULL })
  if (!is.null(ev)) {
    write.table(ev, out("_evanno.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    message("Evanno table: ", out("_evanno.txt"))
  }
}

# ---- align replicates, merge, plot ----
# replicates of the same K can have the clusters in a different order
# ("label switching"): alignK puts them in the same order, then the
# replicates of each K are averaged (done here because pophelper's mergeQ
# fails with recent R versions)
aligned <- alignK(slist)
k_run   <- sapply(aligned, ncol)
merged  <- lapply(sort(unique(k_run)), function(k) {
  runs <- aligned[k_run == k]
  m <- Reduce(`+`, lapply(runs, as.matrix)) / length(runs)
  m <- as.data.frame(m / rowSums(m))
  rownames(m) <- inds$id
  attr(m, "ind") <- nrow(m)          # attributes required by pophelper
  attr(m, "k") <- k
  m
})
names(merged) <- paste0("K", sort(unique(k_run)))
ks <- sapply(merged, ncol)
sel <- which(ks %in% K_PLOT)
if (length(sel) == 0) stop("None of K_PLOT found in the results (K available: ", paste(sort(ks), collapse = ", "), ")")
merged <- alignK(merged[sel], type = "across")

plot_args <- list(qlist = merged, imgoutput = "join", grplab = inds["pop"],
                  grplabangle = 90, grplabsize = 3, showindlab = FALSE,
                  showsp = TRUE, splab = paste0("K=", ks[sel]),
                  exportpath = OUTDIR,
                  outputfilename = paste0(PREFIX, "_K", min(ks[sel]), "-", max(ks[sel])),
                  imgtype = "pdf", width = 25, height = 3)   # cm; height is per K
if (!is.null(COLORS)) plot_args$clustercol <- COLORS
do.call(plotQ, plot_args)

message("Plots written to ", OUTDIR, "/")
