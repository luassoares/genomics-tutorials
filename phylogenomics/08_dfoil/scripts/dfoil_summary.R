# Summarize the Dfoil results (ExDFOIL) by population.
#
# Usage (in the results/ folder, after associate.sh):
#   Rscript ../scripts/dfoil_summary.R all.summary.alt_appended.txt
#
# Needs introg_classifier.R (from ExDFOIL/Scripts) in the same folder.
# Output: freq_introgression.csv (number of tests in each introgression class)
#         and a table on the screen.

args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) > 0) args[1] else "all.summary.alt_appended.txt"

MyDFOIL <- read.table(infile, header = FALSE)
colnames(MyDFOIL) <- c("ind1", "ind2", "ind3", "ind4",
                       "count1", "count2", "count3", "count4",
                       "introg",
                       "DFO_stat", "DFO_p", "DIL_stat", "DIL_p",
                       "DFI_stat", "DFI_p", "DOL_stat", "DOL_p",
                       "name1", "pop1", "species1",
                       "name2", "pop2", "species2",
                       "name3", "pop3", "species3",
                       "name4", "pop4", "species4")

# Translate the Dfoil code (e.g. "123", "31") into a class with the
# population names: column 9 = Dfoil result, 19/22/25/28 = population of
# P1/P2/P3/P4. (Use 20, 23, 26, 29 to classify by species instead.)
source("introg_classifier.R")
MyDFOIL$introg_class <- unlist(apply(MyDFOIL, 1, introg_classifier, 9, 19, 22, 25, 28))

freq <- sort(table(MyDFOIL$introg_class), decreasing = TRUE)
freq_df <- data.frame(class = names(freq),
                      n_tests = as.vector(freq),
                      percent = round(100 * as.vector(freq) / nrow(MyDFOIL), 2))
write.csv(freq_df, "freq_introgression.csv", row.names = FALSE)

cat("Total tests:", nrow(MyDFOIL), "\n")
cat("Tests with a signal of introgression:", sum(MyDFOIL$introg_class != "none"), "\n\n")
print(head(freq_df, 30), row.names = FALSE)
