# conStruct: plots of the results of run_construct.R.
# Run in RStudio, line by line (see the tutorial).

library(conStruct)

setwd("~/my_project/08_construct")      # change to your folder

# 1. Cross-validation ------------------------------------------------------------
sp  <- as.matrix(read.table("results/my_data_sp_xval_results.txt", header = TRUE))    # spatial
nsp <- as.matrix(read.table("results/my_data_nsp_xval_results.txt", header = TRUE))  # non-spatial
# rows = K (1 to 5), columns = replicates

# mean predictive accuracy for each K (higher = better)
plot(rowMeans(sp), pch = 19, col = "blue", ylim = range(sp, nsp),
     xlab = "K", ylab = "predictive accuracy", main = "spatial (blue) and non-spatial (green)")
points(rowMeans(nsp), pch = 19, col = "green")

# 95% confidence intervals of the spatial model
sp_CI <- apply(sp, 1, function(x) mean(x) + c(-1.96, 1.96) * sd(x) / sqrt(length(x)))
plot(rowMeans(sp), pch = 19, col = "blue", ylim = range(sp_CI),
     xlab = "K", ylab = "predictive accuracy", main = "spatial model, 95% CI")
segments(x0 = 1:nrow(sp), y0 = sp_CI[1, ], x1 = 1:nrow(sp), y1 = sp_CI[2, ], col = "blue", lwd = 2)

# Is K = 2 better than K = 1? (paired by replicate)
t.test(sp[2, ], sp[1, ], paired = TRUE, alternative = "greater")
# Is the spatial model better than the non-spatial one at K = 1?
t.test(sp[1, ], nsp[1, ], paired = TRUE, alternative = "greater")

# 2. Layer contributions (spatial runs K = 1 to 5) -------------------------------
layer_contributions <- matrix(0, nrow = 5, ncol = 5)

load("results/my_data_spK1_conStruct.results.Robj")
load("results/my_data_spK1_data.block.Robj")
layer_contributions[1, 1] <- calculate.layer.contribution(conStruct.results[[1]], data.block)
previous <- conStruct.results[[1]]$MAP$admix.proportions

for (k in 2:5) {
  load(paste0("results/my_data_spK", k, "_conStruct.results.Robj"))
  load(paste0("results/my_data_spK", k, "_data.block.Robj"))
  # the layers can be in a different order in each run: match them to the previous K
  layer_order <- match.layers.x.runs(previous, conStruct.results[[1]]$MAP$admix.proportions)
  layer_contributions[1:k, k] <- calculate.layer.contribution(conStruct.results[[1]], data.block,
                                                               layer.order = layer_order)
  previous <- conStruct.results[[1]]$MAP$admix.proportions[, layer_order]
}

colors <- c("blue", "red", "goldenrod1", "forestgreen", "darkorchid1")
barplot(layer_contributions, col = colors, names.arg = paste0("K=", 1:5),
        ylab = "layer contributions")
round(layer_contributions, 3)

# 3. Admixture plots for the chosen K --------------------------------------------
K <- 3                                       # change to your K
load(paste0("results/my_data_spK", K, "_conStruct.results.Robj"))
admix <- conStruct.results[[1]]$MAP$admix.proportions

load(paste0("results/my_data_spK", K, "_data.block.Robj"))
coords <- data.block$coords

make.structure.plot(admix.proportions = admix, layer.colors = colors[1:K])     # bar plot
make.admix.pie.plot(admix.proportions = admix, coords = coords,
                    layer.colors = colors[1:K], radii = 2.5)                  # pies on the map
