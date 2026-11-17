# Make the EEMS input files .diffs and .coord from the VCF.
# Run in RStudio, line by line (see section 3 of the tutorial).
# The .diffs calculation follows the str2diffs page of the EEMS documentation.

library(vcfR)

setwd("~/my_project/09_eems")            # change to your folder

# 1. Genotypes: number of copies of the alternative allele (0, 1, 2) --------
vcf <- read.vcfR("data/my_data-miss0.9-maf0.01-th100.recode.vcf")
gt  <- extract.gt(vcf)                    # genotypes "0/0", "0/1", "1/1"
gt  <- gsub("|", "/", gt, fixed = TRUE)
geno <- matrix(NA, nrow(gt), ncol(gt))
geno[gt == "0/0"] <- 0
geno[gt == "0/1" | gt == "1/0"] <- 1
geno[gt == "1/1"] <- 2
geno <- t(geno)                           # individuals x SNPs
samples <- colnames(gt)

nIndiv <- nrow(geno)
nSites <- ncol(geno)
nIndiv; nSites                            # write these numbers in the params file

# 2. Average genetic differences between all pairs of individuals ------------
# (bed2diffs_v2: missing genotypes are replaced by the mean of the SNP)
bed2diffs_v2 <- function(Geno) {
  nIndiv <- nrow(Geno)
  nSites <- ncol(Geno)
  Miss <- is.na(Geno)
  Mean <- matrix(colMeans(Geno, na.rm = TRUE), nrow = nIndiv, ncol = nSites, byrow = TRUE)
  Mean[!Miss] <- 0
  Geno[Miss] <- 0
  Geno <- Geno + Mean
  Sim <- Geno %*% t(Geno) / nSites
  SelfSim <- diag(Sim)
  vector1s <- rep(1, nIndiv)
  SelfSim %*% t(vector1s) + vector1s %*% t(SelfSim) - 2 * Sim
}

diffs <- round(bed2diffs_v2(geno), 6)

# check: a valid matrix has exactly one positive eigenvalue
sum(eigen(diffs)$values > 0)              # must be 1

write.table(diffs, "data/my_data.diffs", col.names = FALSE, row.names = FALSE, quote = FALSE)

# 3. Coordinates in the same order as the VCF (longitude latitude) ----------
coords <- read.table("data/coords.txt", header = FALSE, col.names = c("id", "lon", "lat"))
coords <- coords[match(samples, coords$id), ]
sum(is.na(coords$lon))                    # must be 0
write.table(coords[, c("lon", "lat")], "data/my_data.coord",
            col.names = FALSE, row.names = FALSE, quote = FALSE)
