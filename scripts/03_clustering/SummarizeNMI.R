###################################################################################
#              Summarize NMI values from random down-sampling analyses            #
###################################################################################

# This Rscript allows to combine NMI values from the random down-sampling analyses and summarize them.
# R code ran on R 4.2.0

# Input files include the normalized mutual information (NMI) score from each random down-sampling analysis.

# Output files include a table with all NMI values as well as a visualization.

# Define variables
Variables <- list(
  
  ## pattern of file to look for
  pattern = "External-validation-ContributionOfEachTumorToEachCluster.txt", 
  
  ## path of directory that contains results from each random down-sampling analysis
  input.path = "."
)

# get all files containing the NMI values
temp <- list.files(path = Variables$input.path, pattern = Variables$pattern)

# get number of sets per gene number
nSets <- gsub(paste0("_", Variables$pattern, ".*"), "", temp)
nSets <- gsub(".*Geneset", "Geneset", nSets)

# get each gene number tested
nGenes <- gsub("Genes_.*", "", temp)
nGenes <- gsub(".*_", "", nGenes)

names(temp) <- paste0(nGenes, "_", nSets)

# create table to store all NMI values
allNMI.tbl <- data.frame(matrix(, nrow=0, ncol=3))
names(allNMI.tbl) <- c("nGenes", "Geneset", "NMI")

for(i in 1:length(temp)){
  
  # Import file containing NMI value
  NMI.value <- read.delim(temp[i], sep = "\t", row.names = NULL, header = FALSE)
  
  # Extract NMI value
  NMI.value <- NMI.value[5, 1]
  NMI.value <- gsub("normalized mutual information  : ", "", NMI.value)
  NMI.value <- as.numeric(gsub(" ", "", NMI.value))
  
  # add to table
  allNMI.tbl[i, "nGenes"] <- as.numeric(gsub("_.*", "", names(temp[i])))
  allNMI.tbl[i, "Geneset"] <- gsub(".*_", "", names(temp[i]))
  allNMI.tbl[i, "NMI"] <- NMI.value
  rm(NMI.value); gc()
}
allNMI.tbl <- allNMI.tbl[order(allNMI.tbl$nGenes, decreasing = FALSE),]
nGenes.order <- unique(allNMI.tbl$nGenes)

# Get summary statistics by gene number
library(dplyr)
summary.results <- allNMI.tbl %>%
  group_by(nGenes) %>%
  summarize(Mean = mean(NMI, na.rm = TRUE), Median = median(NMI, na.rm = TRUE), SD = sd(NMI, na.rm = TRUE))
write.csv(summary.results, paste0(Sys.Date(), "_SummarizedNMIacrossRandomDownsamplingAnalyses.csv"))

# Visualize results
library(ggplot2)
pdf(paste0(Sys.Date(), "_PlotShowingSummarizedNMIacrossRandomDownsamplingAnalyses.pdf"))
ggplot(allNMI.tbl, aes(x = factor(nGenes), y = NMI)) +
  geom_boxplot() +
  theme_classic() +
  theme(legend.position = "none") +
  xlab("")
dev.off()