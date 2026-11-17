# Make a triangular grid for FEEMS that covers only the study area, finer than
# the grids that come with FEEMS. Run in RStudio or with Rscript (see section 3).

library(dggridR)
library(sf)
sf_use_s2(FALSE)

setwd("~/my_project/10_feems")           # change to your folder

# polygon of the study area (longitude latitude), the same as for EEMS
outer <- read.table("data/my_data.outer", col.names = c("lon", "lat"))

# resolution of the grid: 6 = ~110 km between nodes, 7 = ~55 km,
#                         8 = ~28 km, 9 = ~14 km
res <- 8

dggs <- dgconstruct(res = res, precision = 30, projection = "ISEA",
                    aperture = 4, topology = "TRIANGLE")
dggetres(dggs)                            # spacing (km) of each resolution

grid <- dgrectgrid(dggs,
                   minlat = min(outer$lat) - 0.5, maxlat = max(outer$lat) + 0.5,
                   minlon = min(outer$lon) - 0.5, maxlon = max(outer$lon) + 0.5,
                   cellsize = 0.01)
nrow(grid)                                # number of triangles

st_write(grid, "data/grid_res8.shp")
