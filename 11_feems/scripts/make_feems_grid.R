#!/usr/bin/env Rscript
# Make a triangular grid (shapefile) for FEEMS covering only the study area,
# at a finer resolution than the global grids that come with FEEMS
# (grid_100 = res 6, ~110 km between nodes; grid_250 = res 5; grid_500 = res 4).
#
# Usage: edit the settings below, then
#   Rscript make_feems_grid.R
#
# Output: GRID_FILE (.shp, .shx, .dbf, .prj), input of run_feems.py (--grid)

# ---- settings ----
OUTER     <- "feems/my_data.outer"   # habitat polygon: longitude latitude, one vertex per line
GRID_FILE <- "feems/grid_res8.shp"
RES       <- 8          # dggridR resolution: 6 ~110 km, 7 ~55 km, 8 ~28 km, 9 ~14 km between nodes
MARGIN    <- 0.5        # degrees added around the polygon

# ---- packages ----
for (p in c("dggridR", "sf")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Install the R package ", p)
}
suppressPackageStartupMessages({ library(dggridR); library(sf) })
sf_use_s2(FALSE)

outer <- read.table(OUTER, header = FALSE)[, 1:2]
colnames(outer) <- c("lon", "lat")
if (any(abs(outer$lat) > 90)) stop("Latitude > 90: the outer file must be longitude latitude")

# same grid type as the FEEMS grids: ISEA projection, aperture 4, triangles
dggs <- dgconstruct(res = RES, precision = 30, projection = "ISEA", aperture = 4, topology = "TRIANGLE")
info <- dggetres(dggs)
message(sprintf("Resolution %d: cell spacing ~%.0f km, cell area ~%.0f km2",
                RES, info$spacing_km[info$res == RES], info$area_km[info$res == RES]))

grid <- dgrectgrid(dggs,
                   minlat = min(outer$lat) - MARGIN, maxlat = max(outer$lat) + MARGIN,
                   minlon = min(outer$lon) - MARGIN, maxlon = max(outer$lon) + MARGIN,
                   cellsize = 0.01)
dir.create(dirname(GRID_FILE), recursive = TRUE, showWarnings = FALSE)
st_write(grid, GRID_FILE, delete_dsn = file.exists(GRID_FILE), quiet = TRUE)
message(sprintf("%d triangles written to %s", nrow(grid), GRID_FILE))
