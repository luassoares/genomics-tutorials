#!/usr/bin/env Rscript
# Make the EEMS input files <datapath>.diffs and <datapath>.coord from a VCF
# and a table of coordinates. Replaces bed2diffs (no plink or plinkio needed).
# The habitat outline (<datapath>.outer) is made by hand (see the tutorial).
#
# Usage: edit the settings below, then
#   Rscript vcf2diffs.R
#
# Output:
#   <DATAPATH>.diffs   matrix of average pairwise genetic differences (nIndiv x nIndiv)
#   <DATAPATH>.coord   longitude <space> latitude, one line per individual, same order as .diffs
#   <DATAPATH>_samples.txt   sample IDs in the order of .diffs and .coord
#   prints nIndiv and nSites for the params file

# ---- settings ----
VCF      <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
COORDS   <- "coords_individuals.txt"  # <sample ID> <tab> <longitude> <tab> <latitude>, no header
DATAPATH <- "eems/data/my_data"       # output path + prefix (the folder is created)
METHOD   <- "v1"                      # "v1": pairwise-complete SNPs (bed2diffs_v1)
                                      # "v2": missing genotypes imputed with the mean (bed2diffs_v2)

# ---- packages ----
if (!requireNamespace("vcfR", quietly = TRUE)) stop("Install the R package vcfR")
suppressPackageStartupMessages(library(vcfR))
dir.create(dirname(DATAPATH), recursive = TRUE, showWarnings = FALSE)

# ---- genotypes: number of copies of the alternative allele (0, 1, 2) ----
vcf <- read.vcfR(VCF, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT")             # SNPs x samples
gt  <- gsub("|", "/", gt, fixed = TRUE)
geno <- matrix(NA_real_, nrow(gt), ncol(gt), dimnames = dimnames(gt))
geno[gt == "0/0"] <- 0
geno[gt %in% c("0/1", "1/0")] <- 1
geno[gt == "1/1"] <- 2
geno <- t(geno)                                     # samples x SNPs
geno <- geno[, colSums(!is.na(geno)) > 0, drop = FALSE]

n_ind  <- nrow(geno)
n_site <- ncol(geno)

# ---- differences (the two functions of the EEMS str2diffs README) ----
diffs_v1 <- function(Geno) {
  # Dij = mean over SNPs genotyped in both i and j of (z_i - z_j)^2
  obs <- !is.na(Geno) * 1
  G0  <- Geno; G0[is.na(G0)] <- 0
  sq  <- G0^2
  n_ij <- obs %*% t(obs)                            # SNPs genotyped in both
  num  <- sq %*% t(obs) + obs %*% t(sq) - 2 * G0 %*% t(G0)
  D <- num / n_ij
  diag(D) <- 0
  D
}
diffs_v2 <- function(Geno) {
  n <- nrow(Geno); s <- ncol(Geno)
  Miss <- is.na(Geno)
  Mean <- matrix(colMeans(Geno, na.rm = TRUE), nrow = n, ncol = s, byrow = TRUE)
  Mean[!Miss] <- 0
  Geno[Miss] <- 0
  Geno <- Geno + Mean
  Sim <- Geno %*% t(Geno) / s
  SelfSim <- diag(Sim)
  ones <- rep(1, n)
  SelfSim %*% t(ones) + ones %*% t(SelfSim) - 2 * Sim
}

D <- if (METHOD == "v1") diffs_v1(geno) else diffs_v2(geno)
if (any(!is.finite(D))) stop("Some pairs of individuals share no genotyped SNP: use METHOD <- \"v2\" or remove them")
D <- round(D, 6)

# A valid dissimilarity matrix has exactly one positive eigenvalue
ev <- eigen(D, symmetric = TRUE, only.values = TRUE)$values
n_pos <- sum(ev > 1e-8)
if (n_pos != 1) {
  warning(sprintf("The diffs matrix has %d positive eigenvalues (should be 1). ", n_pos),
          if (METHOD == "v1") "Try METHOD <- \"v2\"." else "Check the data.")
}

write.table(D, paste0(DATAPATH, ".diffs"), col.names = FALSE, row.names = FALSE, quote = FALSE)

# ---- coordinates in the same order ----
co <- read.table(COORDS, sep = "\t", header = FALSE, colClasses = c("character", "numeric", "numeric"),
                 col.names = c("id", "lon", "lat"), strip.white = TRUE)
miss <- setdiff(rownames(geno), co$id)
if (length(miss) > 0) stop("Samples without coordinates: ", paste(miss, collapse = ", "))
co <- co[match(rownames(geno), co$id), ]
if (any(abs(co$lat) > 90)) stop("Latitude > 90: are the columns in the order ID, longitude, latitude?")
write.table(co[, c("lon", "lat")], paste0(DATAPATH, ".coord"), col.names = FALSE, row.names = FALSE, quote = FALSE)
writeLines(rownames(geno), paste0(DATAPATH, "_samples.txt"))

message(sprintf("Written %s.diffs and %s.coord", DATAPATH, DATAPATH))
message(sprintf("nIndiv = %d\nnSites = %d", n_ind, n_site))
if (!file.exists(paste0(DATAPATH, ".outer"))) {
  message("Now make ", DATAPATH, ".outer (habitat polygon, see the tutorial).")
}
