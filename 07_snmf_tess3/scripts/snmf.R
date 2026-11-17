# sNMF (LEA): ancestry coefficients, bar plot and map.
# Run in RStudio, line by line (see the tutorial).

library(vcfR)
library(LEA)
library(ggplot2)
library(scatterpie)

setwd("~/my_project/07_snmf_tess3")      # change to your folder

# 1. Read the data -----------------------------------------------------------
vcf <- read.vcfR("data/my_data-miss0.9-maf0.01-th100.recode.vcf")
samples <- colnames(vcf@gt)[-1]           # sample IDs, in the VCF order

popmap <- read.table("data/popmap.txt", header = FALSE, col.names = c("pop", "id"))
coords <- read.table("data/coords.txt", header = FALSE, col.names = c("id", "lon", "lat"))

# put the population and coordinates in the same order as the VCF
pop    <- popmap$pop[match(samples, popmap$id)]
coords <- coords[match(samples, coords$id), ]
sum(is.na(pop)); sum(is.na(coords$lon))   # both must be 0

# 2. Convert the VCF to the LEA .geno format ---------------------------------
# one line per SNP, one number per individual: 0, 1, 2 = copies of the
# alternative allele, 9 = missing
gt <- extract.gt(vcf)                     # genotypes "0/0", "0/1", "1/1"
gt <- gsub("|", "/", gt, fixed = TRUE)
geno <- matrix(9, nrow(gt), ncol(gt))
geno[gt == "0/0"] <- 0
geno[gt == "0/1" | gt == "1/0"] <- 1
geno[gt == "1/1"] <- 2
write.table(geno, "results/my_data.geno", col.names = FALSE, row.names = FALSE, sep = "")

# 3. Run sNMF for K = 1 to 10, 10 repetitions -------------------------------
project <- snmf("results/my_data.geno", K = 1:10, repetitions = 10,
                entropy = TRUE, alpha = 100, CPU = 4, project = "new")
# to open the results again later, without running:
# project <- load.snmfProject("results/my_data.snmfProject")

# 4. Choose K: cross-entropy ------------------------------------------------
plot(project, col = "blue", pch = 19, cex = 1.2)      # lowest cross-entropy = best K

# 5. Ancestry coefficients for the chosen K ---------------------------------
K <- 3                                                # change to your K
best_run <- which.min(cross.entropy(project, K = K))  # run with the lowest cross-entropy
q <- Q(project, K = K, run = best_run)
colnames(q) <- paste0("C", 1:K)

q_table <- data.frame(individual = samples, population = pop,
                      lon = coords$lon, lat = coords$lat, q)
q_table$cluster <- colnames(q)[apply(q, 1, which.max)]    # most likely cluster
write.csv(q_table, "results/snmf_Q_K3.csv", row.names = FALSE)
table(q_table$population, q_table$cluster)

# 6. Bar plot -----------------------------------------------------------------
colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

q_long <- data.frame(individual = rep(samples, K), population = rep(pop, K),
                     cluster = rep(colnames(q), each = length(samples)),
                     q = as.vector(q))
q_long$individual <- factor(q_long$individual, levels = samples[order(pop)])   # group by population

ggplot(q_long, aes(x = individual, y = q, fill = cluster)) +
  geom_col(width = 1) +
  facet_grid(~ population, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        strip.text.x = element_text(angle = 90))
ggsave("results/snmf_barplot_K3.pdf", width = 12, height = 4)

# 7. Map: one pie per site (mean of the individuals of the same site) -------
sites <- aggregate(q_table[, colnames(q)], by = list(lon = q_table$lon, lat = q_table$lat), FUN = mean)
world <- map_data("world")

ggplot() +
  geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey95", colour = "grey70") +
  geom_scatterpie(data = sites, aes(x = lon, y = lat, r = 0.05), cols = colnames(q)) +
  scale_fill_manual(values = colors) +
  coord_quickmap(xlim = range(sites$lon) + c(-0.5, 0.5), ylim = range(sites$lat) + c(-0.5, 0.5)) +
  labs(x = "Longitude", y = "Latitude", fill = "Cluster") +
  theme_bw()
ggsave("results/snmf_map_K3.pdf", width = 8, height = 7)
# r = 0.05 is the size of the pies, in degrees: change it to fit your map
