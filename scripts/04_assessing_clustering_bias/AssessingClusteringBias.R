##########################################################################################################
#              Assessing whether cell clusters are enriched for cells from a specific feature            #
##########################################################################################################

# This Rscript allows to evaluate the contribution of different cell features to the identified cell clusters.
# R code ran on R 4.2.0

# Input files include clustering results and cells' annotations. Clustering results can be obtained from the CellGroupingAnalysis_HCPC_RunCode.R

# Output files include the normalized mutual information (NMI) score, other clustering metrics as well as an alluvial plot as visualization.

# Define variables
Variables <- list(
  
  ## path to clustering results and cell features. Table in csv format with cells as rows and each feature as column
  cluster.path = "2026-07-25_DatasetTest_NoDataManipulation_MetadataWithClusterIDs_.csv", 
  
  ## name of column containing cell groups (from clustering analysis)
  cell.group.col = "Cluster", 
  
  ## name of column containing the cell feature to assess. Cell feature should be a qualitative label
  feature.compare = "Location", # Sample.name, Location
  
  ## name to add to output file
  out.prefix = "Region", # Tumor, Region
  
  # Color to assign to each group. List of comma-separated colors.
  Color = "black,green,orange,pink"
  
)

# import cells’ annotations
metadata <- read.csv(Variables$cluster.path, header=TRUE, row.names=1)

# Visualize clustering bias
## prepare data
library(dplyr)
data <- table(metadata[, Variables$feature.compare], metadata[, Variables$cell.group.col]) %>% as.data.frame()
names(data) <- c(Variables$feature.compare, Variables$cell.group.col, "CellNumber")

## plot
library(ggplot2)
library(ggalluvial)
colorPalette  <- unlist(strsplit(Variables$Color, split = ","))
pdf(paste0(Sys.Date(), "_AlluvialPlot_ContributionOfEach", Variables$out.prefix, "ToEachCluster.pdf"), width = 3, height = 3)
ggplot(data,
       aes(axis1 = data[,Variables$cell.group.col],
           axis2 = data[,Variables$feature.compare],
           y = CellNumber)) +
  geom_alluvium(aes(fill = data[,Variables$feature.compare]), alpha = 0.8, width = 0.5) +
  geom_stratum(width = 0.5) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values = colorPalette) +
  theme(legend.position = "bottom") +
  labs(fill = Variables$feature.compare)
dev.off()

# Quantify clustering bias
## Compute the normalized mutual information (NMI) between cell feature and the assigned cell clusters
library(ClusterR)
sink(paste0(Sys.Date(), "_External-validation-ContributionOfEach", Variables$out.prefix, "ToEachCluster.txt"))

external_validation(
  as.numeric(factor(metadata[, Variables$feature.compare])),
  as.numeric(factor(metadata[, Variables$cell.group.col])),
  summary_stats = TRUE
)

sink()