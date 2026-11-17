#!/usr/bin/env Rscript
# Ancestry coefficients with sNMF (LEA) and TESS3 (tess3r), from a VCF.
#
# Usage: edit the settings below, then run twice:
#   1. K_PLOT <- NA : converts the VCF, runs sNMF and TESS3 for K = 1..K_MAX and
#                     writes the plots to choose K (cross-entropy, cross-validation). Stops.
#   2. K_PLOT set   : loads the saved runs (no re-run) and writes the Q tables,
#                     bar plots and pie maps for that K.
#
#   Rscript snmf_tess3.R
#
# Input:
#   - VCF of biallelic SNPs
#   - POPMAP: <population> <tab> <sample ID>, no header
#   - COORDS: <sample ID> <tab> <longitude> <tab> <latitude>, decimal degrees, no header
#
# Output (in OUTDIR, file names start with PREFIX):
#   .geno                         genotypes in LEA format (0/1/2, 9 = missing)
#   .snmfProject, .snmf/          sNMF runs (LEA project)
#   _tess3.rds                    TESS3 runs
#   _snmf_cross-entropy.pdf/.txt  cross-entropy per K (pass 1)
#   _tess3_cross-validation.pdf   cross-validation score per K (pass 1)
#   _<method>_K<K>_Q.csv          ancestry coefficients, likely cluster and its probability (pass 2)
#   _<method>_K<K>_barplot.pdf    ancestry bar plot, individuals grouped by population (pass 2)
#   _<method>_K<K>_map.pdf        pie charts on the map, one pie per site (pass 2)

# ---- settings ----
VCF      <- "filtered/my_data-miss9-maf1-th100.recode.vcf"
POPMAP   <- "popmap.txt"
COORDS   <- "coords_individuals.txt"
PREFIX   <- "my_data"
OUTDIR   <- "snmf_tess3"

K_MAX    <- 16             # test K = 1..K_MAX
REPS     <- 10             # replicate runs per K
ALPHA    <- 100            # sNMF regularisation (higher = smoother ancestry; 10-100 for small datasets)
N_CORES  <- 4
RUN_SNMF  <- TRUE
RUN_TESS3 <- TRUE
RERUN    <- FALSE          # TRUE = run again even if saved runs exist

K_PLOT   <- NA             # K for the tables and plots (NA = stop after pass 1)
COLORS   <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
              "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78",
              "#98df8a", "#ff9896", "#c5b0d5", "#c49c94")
PIE_SCALE <- 1.5           # size of the pies on the map
MAP_MARGIN <- 0.3          # degrees added around the sites

# ---- packages ----
needed <- c("vcfR", "LEA", "ggplot2", "maps")
if (RUN_TESS3) needed <- c(needed, "tess3r")
if (!is.na(K_PLOT)) needed <- c(needed, "scatterpie")
for (p in needed) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install the R package ", p, " (see 00_setup)")
}
suppressPackageStartupMessages({
  library(vcfR); library(LEA); library(ggplot2)
})
dir.create(OUTDIR, showWarnings = FALSE)
out <- function(x) file.path(OUTDIR, paste0(PREFIX, x))
if (!is.na(K_PLOT) && K_PLOT > length(COLORS)) stop("Give at least K_PLOT colors")

# ---- read data ----
vcf <- read.vcfR(VCF, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT")                      # SNPs x samples
gt  <- gsub("|", "/", gt, fixed = TRUE)
geno <- matrix(9L, nrow(gt), ncol(gt))                     # 9 = missing (LEA format)
geno[gt == "0/0"] <- 0L
geno[gt %in% c("0/1", "1/0")] <- 1L
geno[gt == "1/1"] <- 2L
samples <- colnames(gt)

popmap <- read.table(POPMAP, sep = "\t", header = FALSE, colClasses = "character",
                     col.names = c("pop", "id"), strip.white = TRUE)
pops <- popmap$pop[match(samples, popmap$id)]
if (any(is.na(pops))) stop("Samples missing from the popmap: ", paste(samples[is.na(pops)], collapse = ", "))

co <- read.table(COORDS, sep = "\t", header = FALSE, colClasses = c("character", "numeric", "numeric"),
                 col.names = c("id", "lon", "lat"), strip.white = TRUE)
miss <- setdiff(samples, co$id)
if (length(miss) > 0) stop("Samples without coordinates: ", paste(miss, collapse = ", "))
co <- co[match(samples, co$id), ]
if (any(abs(co$lat) > 90)) stop("Latitude > 90: are the columns in the order ID, longitude, latitude?")
coords <- as.matrix(co[, c("lon", "lat")])
message(sprintf("%d individuals, %d SNPs, %d populations", length(samples), nrow(geno), length(unique(pops))))

# ---- sNMF ----
geno_file <- out(".geno")
proj_file <- out(".snmfProject")
snmf_proj <- NULL
if (RUN_SNMF) {
  if (file.exists(proj_file) && !RERUN) {
    snmf_proj <- load.snmfProject(proj_file)
    message("sNMF: loaded ", proj_file)
  } else {
    # .geno: one line per SNP, one character per individual, no separator
    writeLines(apply(geno, 1, paste, collapse = ""), geno_file)
    message("sNMF: K = 1..", K_MAX, ", ", REPS, " repetitions. Start: ", Sys.time())
    snmf_proj <- snmf(geno_file, K = 1:K_MAX, repetitions = REPS, entropy = TRUE,
                      alpha = ALPHA, CPU = N_CORES, project = "new")
  }
  ce <- t(sapply(1:K_MAX, function(k) {
    x <- cross.entropy(snmf_proj, K = k)
    c(K = k, min = min(x), mean = mean(x), max = max(x))
  }))
  write.table(round(ce, 5), out("_snmf_cross-entropy.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  pdf(out("_snmf_cross-entropy.pdf"), width = 6, height = 4.5)
  plot(ce[, "K"], ce[, "mean"], type = "b", pch = 19, col = "blue", ylim = range(ce[, c("min", "max")]),
       xlab = "Number of ancestral populations (K)", ylab = "Cross-entropy (mean, min-max)")
  segments(ce[, "K"], ce[, "min"], ce[, "K"], ce[, "max"], col = "blue")
  invisible(dev.off())
}

# ---- TESS3 ----
tess_file <- out("_tess3.rds")
tess_obj <- NULL
if (RUN_TESS3) {
  suppressPackageStartupMessages(library(tess3r))
  if (file.exists(tess_file) && !RERUN) {
    tess_obj <- readRDS(tess_file)
    message("TESS3: loaded ", tess_file)
  } else {
    X <- t(geno); X[X == 9L] <- NA                         # individuals x SNPs
    message("TESS3: K = 1..", K_MAX, ", ", REPS, " repetitions. Start: ", Sys.time())
    tess_obj <- tess3(X = X, coord = coords, K = 1:K_MAX, rep = REPS, method = "projected.ls",
                      ploidy = 2, openMP.core.num = N_CORES, keep = "best", mask = 0.05)
    saveRDS(tess_obj, tess_file)
  }
  pdf(out("_tess3_cross-validation.pdf"), width = 6, height = 4.5)
  plot(tess_obj, pch = 19, col = "blue", xlab = "Number of ancestral populations (K)",
       ylab = "Cross-validation score")
  invisible(dev.off())
}

if (is.na(K_PLOT)) {
  message("Choose K from ", OUTDIR, "/*cross-entropy.pdf and *cross-validation.pdf, set K_PLOT and run again.")
  quit(save = "no")
}

# ---- tables and plots for K_PLOT ----
suppressPackageStartupMessages(library(scatterpie))
cl_names <- paste0("C", 1:K_PLOT)
cols <- setNames(COLORS[1:K_PLOT], cl_names)
pop_levels <- unique(popmap$pop[popmap$pop %in% pops])
world <- map_data("world")

plot_method <- function(q, method) {
  colnames(q) <- cl_names
  df <- data.frame(individual = samples, population = pops, lon = co$lon, lat = co$lat,
                   likely_cluster = cl_names[apply(q, 1, which.max)],
                   assignment_prob = round(apply(q, 1, max), 4), round(q, 4))
  write.csv(df, out(sprintf("_%s_K%d_Q.csv", method, K_PLOT)), row.names = FALSE)

  # bar plot: individuals grouped by population, sorted by main cluster within population
  long <- do.call(rbind, lapply(cl_names, function(cn)
    data.frame(individual = df$individual, population = df$population, cluster = cn, q = df[[cn]])))
  ord <- order(factor(df$population, levels = pop_levels), df$likely_cluster, -df$assignment_prob)
  long$individual <- factor(long$individual, levels = df$individual[ord])
  long$population <- factor(long$population, levels = pop_levels)
  p_bar <- ggplot(long, aes(x = individual, y = q, fill = cluster)) +
    geom_col(width = 1) +
    facet_grid(~ population, scales = "free_x", space = "free_x", switch = "x") +
    scale_fill_manual(values = cols) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Ancestry proportion", fill = "Cluster",
         title = sprintf("%s, K = %d", method, K_PLOT)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(),
          panel.spacing.x = unit(0.05, "lines"),
          strip.text.x = element_text(angle = 90, hjust = 1), strip.placement = "outside")
  ggsave(out(sprintf("_%s_K%d_barplot.pdf", method, K_PLOT)), p_bar, width = 12, height = 4)

  # map: one pie per site (mean Q of the individuals at the same coordinates)
  site <- paste(df$lon, df$lat)
  sites <- aggregate(df[, cl_names], by = list(lon = df$lon, lat = df$lat), FUN = mean)
  sites$n <- as.vector(table(site)[paste(sites$lon, sites$lat)])
  xl <- range(sites$lon) + c(-1, 1) * MAP_MARGIN
  yl <- range(sites$lat) + c(-1, 1) * MAP_MARGIN
  r <- PIE_SCALE * 0.02 * max(diff(xl), diff(yl))         # pie radius in degrees
  p_map <- ggplot() +
    geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", colour = "grey70",
                 linewidth = 0.2) +
    geom_scatterpie(data = sites, aes(x = lon, y = lat, r = r), cols = cl_names, colour = "grey30") +
    scale_fill_manual(values = cols) +
    coord_quickmap(xlim = xl, ylim = yl) +
    labs(x = "Longitude", y = "Latitude", fill = "Cluster",
         title = sprintf("%s, K = %d (one pie per site)", method, K_PLOT)) +
    theme_bw(base_size = 11) + theme(panel.grid = element_blank())
  ggsave(out(sprintf("_%s_K%d_map.pdf", method, K_PLOT)), p_map, width = 8, height = 7)

  message(method, ": populations (rows) x most likely cluster (columns)")
  print(table(population = factor(df$population, levels = pop_levels), cluster = df$likely_cluster))
}

if (!is.null(snmf_proj)) {
  best <- which.min(cross.entropy(snmf_proj, K = K_PLOT))   # run with the lowest cross-entropy
  plot_method(LEA::Q(snmf_proj, K = K_PLOT, run = best), "sNMF")
}
if (!is.null(tess_obj)) {
  plot_method(tess3r::qmatrix(tess_obj, K = K_PLOT), "TESS3")
}
message("Results written to ", OUTDIR, "/")
