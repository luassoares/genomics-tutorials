#!/usr/bin/env Rscript
# Detect outlier SNPs (candidates for selection) with pcadapt and write their
# positions, so they can be removed with vcftools to get a VCF without them (optional, see the tutorial).
#
# Usage:
#   1. First run with K_FINAL <- NA: makes the scree plot and the PCA plots
#      (pcadapt_screeplot.svg, pcadapt_scores_PC*.svg) and stops.
#   2. Look at the scree plot, choose K (see the tutorial), set K_FINAL and
#      run again: writes the outliers.
#
#   Rscript pcadapt_outliers.R
#
# Output (in OUTDIR):
#   <prefix>_screeplot.svg        proportion of variance explained by each PC
#   <prefix>_scores_PC1-2.svg ... individuals on the PCs, colored by population
#   <prefix>_manhattan.svg        -log10(p-value) of each SNP
#   <prefix>_pvalues.svg          histogram of the p-values
#   <prefix>_outliers.txt         CHROM <tab> POS of the outlier SNPs (input of vcftools --exclude-positions)

# ---- settings ----
VCF      <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
POPMAP   <- "popmap.txt"   # <population> <tab> <sample ID>, no header
PREFIX   <- "my_data_pcadapt"
OUTDIR   <- "pcadapt"
K_MAX    <- 20             # number of PCs for the scree plot
K_FINAL  <- NA             # number of PCs to use for the outlier test (NA = only make the plots)
ALPHA    <- 0.05           # significance level after correction
CORRECTION <- "bonferroni" # "bonferroni" (strict) or "BH" (false discovery rate)

# ---- packages ----
for (p in c("pcadapt", "vcfR")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install the R package ", p, " (see 00_setup)")
}
suppressPackageStartupMessages({
  library(pcadapt)
  library(vcfR)
})
dir.create(OUTDIR, showWarnings = FALSE)
out <- function(x) file.path(OUTDIR, paste0(PREFIX, x))

# ---- read data ----
# pcadapt reads only the genotypes; vcfR is used to get the sample order and
# the CHROM and POS of each SNP (in the same order as the pcadapt matrix).
vcf  <- read.vcfR(VCF, verbose = FALSE)
fix  <- getFIX(vcf)
samples <- colnames(vcf@gt)[-1]
geno <- read.pcadapt(VCF, type = "vcf")

popmap <- read.table(POPMAP, sep = "\t", header = FALSE, colClasses = "character",
                     col.names = c("pop", "id"), strip.white = TRUE)
pop <- popmap$pop[match(samples, popmap$id)]
if (any(is.na(pop))) {
  stop("Samples missing from the popmap: ", paste(samples[is.na(pop)], collapse = ", "))
}

# ---- step 1: scree plot and PCA ----
x <- pcadapt(input = geno, K = K_MAX)

# the outlier positions are taken from the VCF by row number, so the two
# must have the same SNPs in the same order
if (length(x$pvalues) != nrow(fix)) {
  stop(sprintf("pcadapt read %d SNPs but the VCF has %d. Use a VCF with only biallelic SNPs.",
               length(x$pvalues), nrow(fix)))
}

svg(out("_screeplot.svg"), width = 7, height = 5)
print(plot(x, option = "screeplot", K = K_MAX))
invisible(dev.off())

for (pcs in list(c(1, 2), c(3, 4), c(5, 6))) {
  if (max(pcs) > K_MAX) next
  svg(out(sprintf("_scores_PC%d-%d.svg", pcs[1], pcs[2])), width = 7, height = 6)
  print(plot(x, option = "scores", i = pcs[1], j = pcs[2], pop = pop))
  invisible(dev.off())
}

if (is.na(K_FINAL)) {
  message("Plots written to ", OUTDIR, "/. Choose K from the scree plot, set K_FINAL and run again.")
  quit(save = "no")
}

# ---- step 2: outlier test with K_FINAL PCs ----
x <- pcadapt(input = geno, K = K_FINAL)

svg(out("_manhattan.svg"), width = 9, height = 4)
print(plot(x, option = "manhattan"))
invisible(dev.off())

svg(out("_pvalues.svg"), width = 6, height = 4)
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
invisible(dev.off())

padj <- p.adjust(x$pvalues, method = CORRECTION)
outliers <- which(padj < ALPHA)

write.table(fix[outliers, c("CHROM", "POS"), drop = FALSE], out("_outliers.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

message(sprintf("K = %d, %s correction, alpha = %g: %d outlier SNPs out of %d (%.1f%%)",
                K_FINAL, CORRECTION, ALPHA, length(outliers), nrow(fix),
                100 * length(outliers) / nrow(fix)))
message("Outlier positions: ", out("_outliers.txt"))
