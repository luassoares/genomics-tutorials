# Optional: find outlier SNPs with pcadapt (see section 5 of the tutorial).
# Run in RStudio, line by line.

library(pcadapt)
library(vcfR)

setwd("~/my_project/03_filtering")    # change to your analysis folder

vcf_file <- "results/my_data-miss0.9-maf0.01-th100.recode.vcf"

# 1. Read the data
geno <- read.pcadapt(vcf_file, type = "vcf")
vcf  <- read.vcfR(vcf_file)
popmap <- read.table("data/popmap.txt", header = FALSE, col.names = c("pop", "id"))
pop <- popmap$pop[match(colnames(vcf@gt)[-1], popmap$id)]   # population of each sample, in the VCF order

# 2. Choose K (number of PCs): scree plot and PCA plots
x <- pcadapt(geno, K = 20)
plot(x, option = "screeplot")
plot(x, option = "scores", pop = pop)                 # PC1 x PC2
plot(x, option = "scores", i = 3, j = 4, pop = pop)   # PC3 x PC4

# 3. Run with the chosen K
K <- 5                                                # change to your K
x <- pcadapt(geno, K = K)
plot(x, option = "manhattan")
hist(x$pvalues, breaks = 50, col = "orange", main = "", xlab = "p-values")

# 4. Outliers (Bonferroni correction)
padj <- p.adjust(x$pvalues, method = "bonferroni")
outliers <- which(padj < 0.05)
length(outliers)

# 5. Write the position (CHROM, POS) of the outliers
positions <- getFIX(vcf)[outliers, c("CHROM", "POS")]
write.table(positions, "results/pcadapt_outliers.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
