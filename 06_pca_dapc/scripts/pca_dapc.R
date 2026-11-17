#!/usr/bin/env Rscript
# PCA and DAPC (adegenet) from a VCF and a popmap.
#
# Usage: edit the settings below, then run in three passes:
#   1. K <- NA            : PCA plots + BIC plot to choose the number of clusters (K). Stops.
#   2. K set, N_PCA_DAPC <- NA : a-score plot to choose the number of PCs for DAPC. Stops.
#   3. K and N_PCA_DAPC set   : DAPC scatter plot, membership plot and table.
#
#   Rscript pca_dapc.R
#
# Output (in OUTDIR, file names start with PREFIX):
#   _PCA_variance.svg          % of variance explained by each PC
#   _PCA_PC1-PC2.svg, _PCA_PC1-PC3.svg   individuals on the PCs, colored by population
#   _PCA_scores.csv            PC scores of each individual
#   _DAPC_BIC.svg              BIC for K = 1..MAX_K (pass 1)
#   _DAPC_ascore.svg           a-score optimisation (pass 2)
#   _DAPC_K<K>_scatter.svg     individuals on the discriminant functions (pass 3)
#   _DAPC_K<K>_compoplot.svg   membership probabilities of each individual (pass 3)
#   _DAPC_K<K>_membership.csv  cluster and membership probabilities of each individual (pass 3)

# ---- settings ----
VCF        <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
POPMAP     <- "popmap.txt"   # <population> <tab> <sample ID>, no header
PREFIX     <- "my_data"
OUTDIR     <- "pca_dapc"
N_PC       <- 10             # PCs kept in the PCA

MAX_K      <- 15             # largest K tested by find.clusters
N_PCA_CLUST <- NA            # PCs used by find.clusters (NA = all, as recommended)
K          <- NA             # number of clusters, chosen from the BIC plot (NA = stop after pass 1)
N_PCA_DAPC <- NA             # PCs used by the DAPC, chosen with the a-score (NA = stop after pass 2)
N_DA       <- NA             # discriminant functions kept (NA = K - 1, the maximum)
SEED       <- 1234           # find.clusters uses random starts: fix the seed to repeat the result

# ---- packages ----
for (p in c("vcfR", "adegenet", "ggplot2")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install the R package ", p, " (see 00_setup)")
}
suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(ggplot2)
})
dir.create(OUTDIR, showWarnings = FALSE)
out <- function(x) file.path(OUTDIR, paste0(PREFIX, x))
set.seed(SEED)

# ---- read data ----
vcf <- read.vcfR(VCF, verbose = FALSE)
gl  <- vcfR2genlight(vcf)

popmap <- read.table(POPMAP, sep = "\t", header = FALSE, colClasses = "character",
                     col.names = c("pop", "id"), strip.white = TRUE)
pops <- popmap$pop[match(indNames(gl), popmap$id)]
if (any(is.na(pops))) {
  stop("Samples missing from the popmap: ", paste(indNames(gl)[is.na(pops)], collapse = ", "))
}
pop(gl) <- pops
message(sprintf("%d individuals, %d SNPs, %d populations", nInd(gl), nLoc(gl), nPop(gl)))

# ---- PCA ----
# glPca replaces missing genotypes by the mean allele frequency
pca <- glPca(gl, nf = N_PC)
var_pc <- 100 * pca$eig / sum(pca$eig)

svg(out("_PCA_variance.svg"), width = 6, height = 4)
barplot(var_pc[1:min(20, length(var_pc))], names.arg = seq_len(min(20, length(var_pc))),
        xlab = "PC", ylab = "% of variance", col = "grey60", border = NA)
invisible(dev.off())

scores <- data.frame(individual = indNames(gl), population = pops, pca$scores)
write.csv(scores, out("_PCA_scores.csv"), row.names = FALSE)

plot_pcs <- function(a, b) {
  ggplot(scores, aes(x = .data[[paste0("PC", a)]], y = .data[[paste0("PC", b)]], fill = population)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(shape = 21, size = 3, alpha = 0.85) +
    labs(x = sprintf("PC%d (%.1f%%)", a, var_pc[a]), y = sprintf("PC%d (%.1f%%)", b, var_pc[b])) +
    theme_bw(base_size = 12)
}
ggsave(out("_PCA_PC1-PC2.svg"), plot_pcs(1, 2), width = 7, height = 5.5)
if (N_PC >= 3) ggsave(out("_PCA_PC1-PC3.svg"), plot_pcs(1, 3), width = 7, height = 5.5)
message("PCA plots written to ", OUTDIR, "/")

# ---- DAPC pass 1: choose K ----
n_pca_clust <- if (is.na(N_PCA_CLUST)) nInd(gl) - 1 else N_PCA_CLUST
clust_all <- find.clusters(gl, n.pca = n_pca_clust, max.n.clust = MAX_K,
                           choose.n.clust = FALSE, criterion = "diffNgroup")
bic <- data.frame(K = seq_along(clust_all$Kstat), BIC = as.numeric(clust_all$Kstat))

svg(out("_DAPC_BIC.svg"), width = 6, height = 4)
plot(bic$K, bic$BIC, type = "b", pch = 19, xlab = "Number of clusters (K)", ylab = "BIC",
     main = "Choose K at the elbow (or the lowest BIC)")
invisible(dev.off())
print(bic, row.names = FALSE)

if (is.na(K)) {
  message("BIC plot: ", out("_DAPC_BIC.svg"), ". Choose K, set K and run again.")
  quit(save = "no")
}

grp <- find.clusters(gl, n.pca = n_pca_clust, n.clust = K)
n_da <- if (is.na(N_DA)) K - 1 else min(N_DA, K - 1)

# ---- DAPC pass 2: choose the number of PCs (a-score) ----
if (is.na(N_PCA_DAPC)) {
  dapc_tmp <- dapc(gl, grp$grp, n.pca = min(n_pca_clust, nInd(gl) - 1), n.da = n_da)
  svg(out("_DAPC_ascore.svg"), width = 7, height = 5)
  asc <- optim.a.score(dapc_tmp)
  invisible(dev.off())
  message("a-score: best number of PCs = ", asc$best,
          ". Set N_PCA_DAPC (e.g. ", asc$best, ") and run again.")
  quit(save = "no")
}

# ---- DAPC pass 3: final DAPC ----
dapc1 <- dapc(gl, grp$grp, n.pca = N_PCA_DAPC, n.da = n_da)
cols <- funky(K)

svg(out(sprintf("_DAPC_K%d_scatter.svg", K)), width = 8, height = 7)
if (n_da >= 2) {
  scatter(dapc1, scree.da = TRUE, bg = "white", cell = 0, cstar = 0, col = cols, solid = 0.8,
          cex = 2, clab = 0, leg = TRUE, txt.leg = paste("Cluster", 1:K),
          posi.da = "bottomright", scree.pca = TRUE, posi.pca = "bottomleft")
} else {
  scatter(dapc1, col = cols, legend = TRUE, txt.leg = paste("Cluster", 1:K))   # K = 2: one function
}
invisible(dev.off())

# membership plot, individuals ordered by population
ord <- order(pops)
svg(out(sprintf("_DAPC_K%d_compoplot.svg", K)), width = 10, height = 4)
par(mar = c(7, 4, 1, 1))
barplot(t(dapc1$posterior[ord, , drop = FALSE]), col = cols, border = NA, space = 0,
        names.arg = pops[ord], las = 2, cex.names = 0.5, ylab = "Membership probability")
invisible(dev.off())

memb <- data.frame(individual = indNames(gl), population = pops,
                   cluster = dapc1$assign, round(dapc1$posterior, 4))
write.csv(memb, out(sprintf("_DAPC_K%d_membership.csv", K)), row.names = FALSE)

message("Populations (rows) x DAPC clusters (columns):")
print(table(population = pops, cluster = dapc1$assign))
message("DAPC results written to ", OUTDIR, "/")
