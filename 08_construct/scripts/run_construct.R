#!/usr/bin/env Rscript
# conStruct: spatial vs non-spatial models of population structure.
#   1. cross-validation for K = 1..K_MAX (spatial and non-spatial models)
#   2. full runs of the spatial (and non-spatial) model for each K, used for
#      the layer contributions and the plots (plot_construct.R)
#
# Usage: edit the settings below, then
#   Rscript run_construct.R          (inside byobu/tmux: it takes hours)
#
# Input:
#   - VCF of biallelic SNPs
#   - COORDS: <sample ID> <tab> <longitude> <tab> <latitude>, decimal degrees, no header
#
# Output (in OUTDIR):
#   <prefix>_freqs.rds, <prefix>_coords.rds, <prefix>_geodist.rds   data used in the runs
#   xval/<prefix>_sp_xval_results.txt, xval/<prefix>_nsp_xval_results.txt   cross-validation
#   runs/<prefix>_spK<K>_conStruct.results.Robj, ..._data.block.Robj    full runs (spatial)
#   runs/<prefix>_nspK<K>_...                                             full runs (non-spatial)

# ---- settings ----
VCF       <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
COORDS    <- "coords_individuals.txt"
PREFIX    <- "my_data"
OUTDIR    <- "construct"
EXCLUDE   <- c()            # sample IDs to remove, e.g. with too much missing data: c("3005924", "3005865")
MAX_MISSING_IND <- 0.5      # remove samples with more than this proportion of missing SNPs

K_MAX     <- 5              # test K = 1..K_MAX
RUN_XVAL  <- TRUE           # step 1: cross-validation
XVAL_REPS <- 10             # cross-validation replicates
XVAL_ITER <- 1e4            # MCMC iterations per cross-validation run
RUN_FULL  <- TRUE           # step 2: full runs for K = 1..K_MAX
FULL_ITER <- 5000           # MCMC iterations of the full runs (increase if the chains do not mix)
N_CORES   <- 16
SEED      <- 1234

# ---- packages ----
for (p in c("conStruct", "vcfR", "geodist", "doParallel", "foreach")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install the R package ", p, " (see 00_setup)")
}
suppressPackageStartupMessages({
  library(conStruct)
  library(vcfR)
  library(geodist)
  library(doParallel)
  library(foreach)
})
set.seed(SEED)
dir.create(file.path(OUTDIR, "xval"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTDIR, "runs"), recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(OUTDIR, paste0(PREFIX, x))

# ---- allele frequencies from the VCF ----
# conStruct needs a matrix samples x loci with the frequency of one allele in
# each sample: 0, 0.5 or 1 for a diploid individual, NA if missing.
vcf <- read.vcfR(VCF, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT")                 # loci x samples
gt  <- gsub("|", "/", gt, fixed = TRUE)
dose <- matrix(NA_real_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
dose[gt == "0/0"] <- 0
dose[gt %in% c("0/1", "1/0")] <- 0.5
dose[gt == "1/1"] <- 1
freqs <- t(dose)                                       # samples x loci

freqs <- freqs[!rownames(freqs) %in% EXCLUDE, , drop = FALSE]
miss_ind <- rowMeans(is.na(freqs))
if (any(miss_ind > MAX_MISSING_IND)) {
  message("Removed (missing data > ", MAX_MISSING_IND, "): ",
          paste(rownames(freqs)[miss_ind > MAX_MISSING_IND], collapse = ", "))
  freqs <- freqs[miss_ind <= MAX_MISSING_IND, , drop = FALSE]
}
# remove loci that became invariant or are missing in all samples
keep <- apply(freqs, 2, function(x) { x <- x[!is.na(x)]; length(x) > 0 && length(unique(x)) > 1 })
freqs <- freqs[, keep, drop = FALSE]

# ---- coordinates and geographic distances ----
co <- read.table(COORDS, sep = "\t", header = FALSE, colClasses = c("character", "numeric", "numeric"),
                 col.names = c("id", "lon", "lat"), strip.white = TRUE)
miss_co <- setdiff(rownames(freqs), co$id)
if (length(miss_co) > 0) stop("Samples without coordinates: ", paste(miss_co, collapse = ", "))
coords <- as.matrix(co[match(rownames(freqs), co$id), c("lon", "lat")])
rownames(coords) <- rownames(freqs)
if (any(abs(coords[, "lat"]) > 90)) stop("Latitude > 90: are the columns in the order ID, longitude, latitude?")

geoDist <- geodist(coords, measure = "geodesic") / 1000   # km

saveRDS(freqs, out("_freqs.rds"))
saveRDS(coords, out("_coords.rds"))
saveRDS(geoDist, out("_geodist.rds"))
message(sprintf("%d samples, %d loci", nrow(freqs), ncol(freqs)))

# ---- parallel backend ----
cl <- makeCluster(N_CORES, type = if (.Platform$OS.type == "unix") "FORK" else "PSOCK")
registerDoParallel(cl)
on.exit(stopCluster(cl), add = TRUE)

# ---- step 1: cross-validation ----
if (RUN_XVAL) {
  message("Cross-validation: K = 1..", K_MAX, ", ", XVAL_REPS, " replicates. Start: ", Sys.time())
  old <- setwd(file.path(OUTDIR, "xval"))
  x.validation(train.prop = 0.9, n.reps = XVAL_REPS, K = 1:K_MAX,
               freqs = freqs, data.partitions = NULL,
               geoDist = geoDist, coords = coords,
               prefix = PREFIX, n.iter = XVAL_ITER,
               make.figs = TRUE, save.files = FALSE,
               parallel = TRUE, n.nodes = N_CORES)
  setwd(old)
  message("Cross-validation done: ", Sys.time())
}

# ---- step 2: full runs ----
if (RUN_FULL) {
  # full path in the prefix: the parallel workers do not follow setwd()
  runs_dir <- normalizePath(file.path(OUTDIR, "runs"))
  runs <- expand.grid(K = 1:K_MAX, spatial = c(TRUE, FALSE))
  foreach(i = seq_len(nrow(runs)), .packages = "conStruct") %dopar% {
    k <- runs$K[i]
    sp <- runs$spatial[i]
    conStruct(spatial = sp, K = k, freqs = freqs, geoDist = geoDist, coords = coords,
              prefix = file.path(runs_dir, sprintf("%s_%sK%d", PREFIX, if (sp) "sp" else "nsp", k)),
              n.chains = 1, n.iter = FULL_ITER, make.figs = TRUE, save.files = TRUE)
    NULL
  }
  message("Full runs done: ", Sys.time())
}
