# Plot the EEMS results with reemsplots2.
# Run in RStudio, line by line (see section 6 of the tutorial).

library(reemsplots2)
library(ggplot2)

setwd("~/my_project/09_eems")            # change to your folder

# 1. Read the results of the three chains --------------------------------------
mcmcpath <- c("results/my_data-chain1",
              "results/my_data-chain2",
              "results/my_data-chain3")
plots <- make_eems_plots(mcmcpath, longlat = TRUE)
names(plots)

# 2. Check the runs ------------------------------------------------------------
plots$pilogl01      # log posterior of each chain: should be flat and similar
plots$rdist01       # observed vs fitted dissimilarities between demes
plots$rdist02       # observed vs fitted dissimilarities within demes
plots$rdist03       # dissimilarity vs geographic distance

# 3. Maps with the sampling sites and country borders ------------------------
coord <- read.table("data/my_data.coord", col.names = c("long", "lat"))
world <- map_data("world")

map_limits <- coord_quickmap(xlim = range(coord$long) + c(-1, 1),
                             ylim = range(coord$lat) + c(-1, 1))

plots$mrates01 +                         # effective migration
  geom_path(data = world, aes(x = long, y = lat, group = group), colour = "grey40") +
  geom_point(data = coord, aes(x = long, y = lat), shape = 1) +
  map_limits
ggsave("results/eems_migration.pdf", width = 7, height = 6)

plots$mrates02 +                         # probability of lower / higher migration
  geom_path(data = world, aes(x = long, y = lat, group = group), colour = "grey40") +
  geom_point(data = coord, aes(x = long, y = lat), shape = 1) +
  map_limits
ggsave("results/eems_migration_probability.pdf", width = 7, height = 6)

plots$qrates01 +                         # effective diversity
  geom_path(data = world, aes(x = long, y = lat, group = group), colour = "grey40") +
  geom_point(data = coord, aes(x = long, y = lat), shape = 1) +
  map_limits
ggsave("results/eems_diversity.pdf", width = 7, height = 6)
