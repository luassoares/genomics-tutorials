# TESS3 (tess3r): ancestry coefficients using the coordinates, bar plot and map.
# Run in RStudio, line by line, after snmf.R (it uses the same .geno file).

library(vcfR)
library(LEA)
library(tess3r)
library(ggplot2)
library(scatterpie)

setwd("~/my_project/07_snmf_tess3")      # change to your folder

# 1. Read the data -----------------------------------------------------------
# genotypes made by snmf.R (9 = missing)
geno <- read.geno("results/my_data.geno")   # individuals x SNPs
geno[geno == 9] <- NA

# sample IDs, populations and coordinates in the VCF order (as in snmf.R)
samples <- colnames(read.vcfR("data/my_data-miss0.9-maf0.01-th100.recode.vcf")@gt)[-1]
popmap <- read.table("data/popmap.txt", header = FALSE, col.names = c("pop", "id"))
coords <- read.table("data/coords.txt", header = FALSE, col.names = c("id", "lon", "lat"))
pop    <- popmap$pop[match(samples, popmap$id)]
coords <- coords[match(samples, coords$id), ]
coord_matrix <- as.matrix(coords[, c("lon", "lat")])

# 2. Run TESS3 for K = 1 to 10 -------------------------------------------------
tess3_obj <- tess3(X = geno, coord = coord_matrix, K = 1:10, rep = 10,
                   method = "projected.ls", ploidy = 2, openMP.core.num = 4,
                   keep = "best", mask = 0.05)
save(tess3_obj, file = "results/tess3_K1-10.Robj")
# to open the results again later: load("results/tess3_K1-10.Robj")

# 3. Choose K: cross-validation score ----------------------------------------
plot(tess3_obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations (K)", ylab = "Cross-validation score")

# 4. Ancestry coefficients for the chosen K ---------------------------------
K <- 3                                                # change to your K
q <- qmatrix(tess3_obj, K = K)
q <- matrix(q, ncol = K)
colnames(q) <- paste0("C", 1:K)

q_table <- data.frame(individual = samples, population = pop,
                      lon = coords$lon, lat = coords$lat, q)
q_table$cluster <- colnames(q)[apply(q, 1, which.max)]
write.csv(q_table, "results/tess3_Q_K3.csv", row.names = FALSE)
table(q_table$population, q_table$cluster)

# 5. Bar plot and map: the same code as in snmf.R (sections 6 and 7) -------
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

q_long <- data.frame(individual = rep(samples, K), population = rep(pop, K),
                     cluster = rep(colnames(q), each = length(samples)),
                     q = as.vector(q))
q_long$individual <- factor(q_long$individual, levels = samples[order(pop)])

ggplot(q_long, aes(x = individual, y = q, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~ population, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_text(angle = 90))
ggsave("results/tess3_barplot_K3.pdf", width = 12, height = 4)

sites <- aggregate(q_table[, colnames(q)], by = list(lon = q_table$lon, lat = q_table$lat), FUN = mean)
world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", colour = "grey70") +
  geom_scatterpie(data = sites, aes(x = lon, y = lat, r = 0.05), cols = colnames(q)) +
  scale_fill_manual(values = colors) +
  coord_quickmap(xlim = range(sites$lon) + c(-0.5, 0.5), ylim = range(sites$lat) + c(-0.5, 0.5)) +
  labs(x = "Longitude", y = "Latitude", fill = "Cluster") +
  theme_bw()
ggsave("results/tess3_map_K3.pdf", width = 8, height = 7)
