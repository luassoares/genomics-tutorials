#!/usr/bin/env Rscript
# Plot the results of run_construct.R:
#   - cross-validation: predictive accuracy of the spatial and non-spatial
#     models for each K, with 95% confidence intervals, and paired t-tests
#   - layer contributions of the spatial runs, K = 1..K_MAX
#   - admixture bar plot and map of pie charts for one K
#
# Usage: edit the settings below, then
#   Rscript plot_construct.R
#
# Output (in OUTDIR/plots, file names start with PREFIX):
#   _xval.pdf                cross-validation, both models
#   _xval_ttests.txt         paired t-tests between consecutive values of K (spatial model)
#   _layer_contributions.pdf / .txt
#   _spK<K>_structure.pdf    admixture proportions per individual
#   _spK<K>_pies.pdf         admixture proportions on the map

# ---- settings ----
OUTDIR  <- "construct"       # same as in run_construct.R
PREFIX  <- "my_data"
K_MAX   <- 5
K_PLOT  <- 3                 # K for the admixture plots
POPMAP  <- "popmap.txt"      # <population> <tab> <sample ID>; used to order the bar plot (NULL = no order)
COLORS  <- c("#3b6a87", "#e07342", "#95c073", "#f8961e", "#8a312b",
             "#735d78", "#5c9ead", "#da8e95", "#7cb684", "#ffe3b2")

# ---- packages ----
if (!requireNamespace("conStruct", quietly = TRUE)) stop("Install the R package conStruct")
suppressPackageStartupMessages(library(conStruct))
pdir <- file.path(OUTDIR, "plots")
dir.create(pdir, showWarnings = FALSE)
out <- function(x) file.path(pdir, paste0(PREFIX, x))
coords <- readRDS(file.path(OUTDIR, paste0(PREFIX, "_coords.rds")))

# ---- cross-validation ----
xfile <- function(m) file.path(OUTDIR, "xval", sprintf("%s_%s_xval_results.txt", PREFIX, m))
if (file.exists(xfile("sp"))) {
  sp  <- as.matrix(read.table(xfile("sp"),  header = TRUE, stringsAsFactors = FALSE))
  nsp <- as.matrix(read.table(xfile("nsp"), header = TRUE, stringsAsFactors = FALSE))
  # rows = K, columns = replicates
  ci <- function(m) apply(m, 1, function(x) mean(x) + c(-1.96, 1.96) * sd(x) / sqrt(length(x)))
  sp_ci <- ci(sp); nsp_ci <- ci(nsp)
  ks <- seq_len(nrow(sp))

  pdf(out("_xval.pdf"), width = 10, height = 4)
  par(mfrow = c(1, 3))
  plot(ks, rowMeans(sp), pch = 19, col = "blue", ylim = range(sp, nsp),
       xlab = "K", ylab = "predictive accuracy", main = "Spatial (blue) vs non-spatial (green)")
  points(ks, rowMeans(nsp), pch = 19, col = "green3")
  for (d in list(list(sp, sp_ci, "blue", "Spatial"), list(nsp, nsp_ci, "green3", "Non-spatial"))) {
    plot(ks, rowMeans(d[[1]]), pch = 19, col = d[[3]], ylim = range(d[[2]]),
         xlab = "K", ylab = "predictive accuracy", main = paste(d[[4]], "(95% CI)"))
    segments(ks, d[[2]][1, ], ks, d[[2]][2, ], col = d[[3]], lwd = 2)
  }
  invisible(dev.off())

  # Is K+1 better than K? Paired by replicate (the same data partition)
  sink(out("_xval_ttests.txt"))
  for (k in ks[-1]) {
    tt <- t.test(sp[k, ], sp[k - 1, ], paired = TRUE, alternative = "greater")
    cat(sprintf("Spatial K=%d vs K=%d: mean difference = %.2f, t = %.2f, p = %.3g\n",
                k, k - 1, tt$estimate, tt$statistic, tt$p.value))
  }
  for (k in ks) {
    tt <- t.test(sp[k, ], nsp[k, ], paired = TRUE, alternative = "greater")
    cat(sprintf("K=%d spatial vs non-spatial: mean difference = %.2f, p = %.3g\n",
                k, tt$estimate, tt$p.value))
  }
  sink()
  cat(readLines(out("_xval_ttests.txt")), sep = "\n")
} else {
  message("No cross-validation results in ", file.path(OUTDIR, "xval"))
}

# ---- layer contributions (spatial runs) ----
rfile <- function(k, what) file.path(OUTDIR, "runs", sprintf("%s_spK%d_%s.Robj", PREFIX, k, what))
layer_contr <- matrix(0, nrow = K_MAX, ncol = K_MAX, dimnames = list(paste0("layer", 1:K_MAX), paste0("K=", 1:K_MAX)))
prev <- NULL
for (k in 1:K_MAX) {
  if (!file.exists(rfile(k, "conStruct.results"))) { message("Missing run K=", k); next }
  load(rfile(k, "conStruct.results"))      # conStruct.results
  load(rfile(k, "data.block"))             # data.block
  res <- conStruct.results[[1]]
  # layers can be in a different order in each run: match them to the previous K
  ord <- if (is.null(prev) || k == 1) seq_len(k) else match.layers.x.runs(prev, res$MAP$admix.proportions)
  layer_contr[1:k, k] <- calculate.layer.contribution(conStruct.results = res, data.block = data.block,
                                                      layer.order = ord)
  prev <- res$MAP$admix.proportions[, ord, drop = FALSE]
}
write.table(round(layer_contr, 4), out("_layer_contributions.txt"), sep = "\t", quote = FALSE)
pdf(out("_layer_contributions.pdf"), width = 6, height = 5)
barplot(layer_contr, col = COLORS[1:K_MAX], ylab = "layer contribution", border = NA)
invisible(dev.off())
print(round(layer_contr, 3))

# ---- admixture plots for K_PLOT ----
if (file.exists(rfile(K_PLOT, "conStruct.results"))) {
  load(rfile(K_PLOT, "conStruct.results"))
  adm <- conStruct.results[[1]]$MAP$admix.proportions
  ids <- rownames(coords)
  ord <- seq_len(nrow(adm))
  if (!is.null(POPMAP)) {
    pm <- read.table(POPMAP, sep = "\t", header = FALSE, colClasses = "character", col.names = c("pop", "id"))
    pops <- pm$pop[match(ids, pm$id)]
    ord <- order(factor(pops, levels = unique(pm$pop)))
  }
  pdf(out(sprintf("_spK%d_structure.pdf", K_PLOT)), width = 12, height = 4)
  make.structure.plot(admix.proportions = adm, sample.order = ord, sample.names = ids,
                      mar = c(6, 4, 2, 2), layer.colors = COLORS[1:K_PLOT])
  invisible(dev.off())
  pdf(out(sprintf("_spK%d_pies.pdf", K_PLOT)), width = 7, height = 7)
  make.admix.pie.plot(admix.proportions = adm, coords = coords, layer.colors = COLORS[1:K_PLOT],
                      radii = 2.5)
  invisible(dev.off())
}
message("Plots written to ", pdir, "/")
