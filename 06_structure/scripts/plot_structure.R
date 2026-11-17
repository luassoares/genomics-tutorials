# Plot STRUCTURE or fastSTRUCTURE results with pophelper.
# Run in RStudio, line by line (see section 6 of the tutorial).

library(pophelper)

setwd("~/my_project/06_structure")      # change to your folder

# 1. Read the result files ---------------------------------------------------
# STRUCTURE: files ending in "_f"; fastSTRUCTURE: files ending in ".meanQ"
files <- list.files("results/structure_out", pattern = "_f$", full.names = TRUE, recursive = TRUE)
# files <- list.files("results/faststructure_out", pattern = "meanQ$", full.names = TRUE, recursive = TRUE)
length(files)

qlist <- readQ(files)

# 2. Summary and Evanno method (STRUCTURE only) ------------------------------
tr <- tabulateQ(qlist)
sr <- summariseQ(tr)
sr                                      # runs and mean likelihood for each K

evanno <- evannoMethodStructure(sr, exportplot = TRUE, exportpath = "results", imgtype = "pdf")
evanno                                  # table: the K with the highest deltaK is the Evanno K

# 3. Population labels -------------------------------------------------------
# one column with header "population", in the same order as the individuals
# in the input file (made by vcf2structure.py)
labels <- read.table("data/my_data_labels.txt", header = TRUE, stringsAsFactors = FALSE)

# 4. Align the runs ----------------------------------------------------------
# the cluster colors can be swapped between runs: align them
qlist_aligned <- alignK(qlist)                    # align the runs of the same K
qlist_aligned <- alignK(qlist_aligned, type = "across")   # same colors for the same cluster across K

# 5. Plot --------------------------------------------------------------------
# which runs to plot: one run per K, K = 2 to 5
ks <- sapply(qlist_aligned, ncol)
sel <- match(2:5, ks)                   # first run of each K

plotQ(qlist_aligned[sel], imgoutput = "join",
      grplab = labels, grplabangle = 90, grplabsize = 3,
      showsp = TRUE, splab = paste0("K=", 2:5),
      showindlab = FALSE,
      exportpath = "results", outputfilename = "structure_K2-5", imgtype = "pdf",
      width = 25, height = 3)

# own colors (one per cluster), e.g.:
# clustercol = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
