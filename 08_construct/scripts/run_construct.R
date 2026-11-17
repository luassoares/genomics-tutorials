# conStruct: cross-validation and runs for K = 1 to 5.
# This takes hours: run it on the server, inside byobu/tmux:
#   conda activate r-popgen
#   Rscript scripts/run_construct.R
# (or open it in RStudio and run it block by block)

library(conStruct)
library(vcfR)
library(geodist)
library(doParallel)

setwd("~/my_project/08_construct")      # change to your folder

# 1. Allele frequencies of each individual from the VCF ---------------------
# conStruct needs a matrix: one row per individual, one column per SNP,
# with the frequency of one allele in the individual (0, 0.5 or 1; NA = missing)
vcf <- read.vcfR("data/my_data-miss0.9-maf0.01-th100.recode.vcf")
gt  <- extract.gt(vcf)                       # genotypes "0/0", "0/1", "1/1"
gt  <- gsub("|", "/", gt, fixed = TRUE)
freqs <- matrix(NA, nrow(gt), ncol(gt))
freqs[gt == "0/0"] <- 0
freqs[gt == "0/1" | gt == "1/0"] <- 0.5
freqs[gt == "1/1"] <- 1
freqs <- t(freqs)                            # individuals x SNPs
rownames(freqs) <- colnames(gt)
dim(freqs)

# Alternative (as in the original analysis): read the STRUCTURE file made in
# 06_structure (by PGDSpider or vcf2structure.py) instead of the VCF.
# The individuals are then in the order of the STRUCTURE file: put the
# coordinates in that order (section 2), using the IDs in its first column.
# freqs <- structure2conStruct(infile = "data/my_data.str", onerowperind = FALSE,
#                              start.loci = 3, start.samples = 2, missing.datum = -9,
#                              outfile = "results/my_data_construct")
# str_ids <- read.table("data/my_data.str", skip = 1)[c(TRUE, FALSE), 1]
# rownames(freqs) <- str_ids

# 2. Coordinates and geographic distances -----------------------------------
coords_table <- read.table("data/coords.txt", header = FALSE, col.names = c("id", "lon", "lat"))
coords_table <- coords_table[match(rownames(freqs), coords_table$id), ]   # same order as the VCF
sum(is.na(coords_table$lon))                 # must be 0
coords  <- as.matrix(coords_table[, c("lon", "lat")])
geoDist <- geodist(coords, measure = "geodesic") / 1000    # distances in km

# 3. Cross-validation: spatial and non-spatial models, K = 1 to 5 -----------
cl <- makeCluster(16, type = "FORK")         # number of processors
registerDoParallel(cl)

xval <- x.validation(train.prop = 0.9,      # 90% of the SNPs to fit, 10% to test
                     n.reps = 10,            # 10 replicates
                     K = 1:5,
                     freqs = freqs,
                     data.partitions = NULL,
                     geoDist = geoDist,
                     coords = coords,
                     prefix = "results/my_data",
                     n.iter = 1e4,
                     make.figs = FALSE,         # TRUE writes plots of every replicate
                     save.files = FALSE,
                     parallel = TRUE,
                     n.nodes = 16)

# 4. Full runs of the spatial model, K = 1 to 5 (for the plots) -------------
foreach(k = 1:5, .packages = "conStruct") %dopar% {
  conStruct(spatial = TRUE, K = k, freqs = freqs, geoDist = geoDist, coords = coords,
            prefix = paste0(getwd(), "/results/my_data_spK", k), n.iter = 5000)
}

stopCluster(cl)
