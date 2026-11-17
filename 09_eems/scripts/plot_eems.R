#!/usr/bin/env Rscript
# Plot the EEMS results of several chains with reemsplots2, adding the
# sampling locations and country borders.
#
# Usage: edit the settings below, then
#   Rscript plot_eems.R
#
# Output (in OUTDIR, file names start with PREFIX):
#   _mrates01.pdf   effective migration surface (log10 m, relative to the mean)
#   _mrates02.pdf   posterior probability that m is above/below the mean
#   _qrates01.pdf   effective diversity surface
#   _qrates02.pdf   posterior probability that q is above/below the mean
#   _rdist01-03.pdf observed vs fitted dissimilarities, and dissimilarity vs distance
#   _pilogl01.pdf   log posterior of each chain (convergence)

# ---- settings ----
MCMCPATH <- Sys.glob("eems/runs/my_data-nDemes300-chain*")   # folders of the chains to combine
COORD    <- "eems/data/my_data.coord"                        # longitude latitude
PREFIX   <- "my_data_eems"
OUTDIR   <- "eems/plots"

# ---- packages ----
for (p in c("reemsplots2", "ggplot2", "maps")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("Install the R package ", p,
         if (p == "reemsplots2") ": remotes::install_github('dipetkov/reemsplots2')" else "")
  }
}
suppressPackageStartupMessages({
  library(reemsplots2)
  library(ggplot2)
})
if (length(MCMCPATH) == 0) stop("No EEMS output folders found")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
out <- function(x) file.path(OUTDIR, paste0(PREFIX, x))

plots <- make_eems_plots(MCMCPATH, longlat = TRUE)

coord <- read.table(COORD, header = FALSE, col.names = c("long", "lat"))
borders <- map_data("world")
lims_x <- range(coord$long) + c(-2, 2)
lims_y <- range(coord$lat) + c(-2, 2)

add_map <- function(p) {
  p +
    geom_path(data = borders, aes(x = long, y = lat, group = group), inherit.aes = FALSE,
              colour = "grey40", linewidth = 0.3) +
    geom_point(data = coord, aes(x = long, y = lat), inherit.aes = FALSE, shape = 1, size = 1.5) +
    coord_quickmap(xlim = lims_x, ylim = lims_y)
}

for (nm in c("mrates01", "mrates02", "qrates01", "qrates02")) {
  # warnings about removed rows are only the parts of the map outside the plotted area
  suppressWarnings(ggsave(out(paste0("_", nm, ".pdf")), add_map(plots[[nm]]), width = 7, height = 6))
}
for (nm in c("rdist01", "rdist02", "rdist03", "pilogl01")) {
  if (!is.null(plots[[nm]])) ggsave(out(paste0("_", nm, ".pdf")), plots[[nm]], width = 6, height = 5)
}
message("Plots written to ", OUTDIR, "/")
