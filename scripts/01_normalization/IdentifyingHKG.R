###################################################################
#          Identification of housekeeping genes (HKG)             #
###################################################################

# This Rscript allows to identify a group of genes that we designated as housekeeping genes (HKG).
    ## We defined HKG according to the following criteria: genes detected in most cells, with relatively constant expression across cells and whose expressions auto-correlate. 
# Input files include the expression data in log2(CPM+1) and the cells' metadata file.
# Output files include the gene clustering results and the identity of the HKG defined.

## Rscript ran on 4.2.0

# Example of command (that generates output files) :
# Rscript ../../IdentifyingHKG.R --file NormalizedExpressionData_OutputStep3.tsv --dataset DatasetTest --analysis_name HKGidentification --corrMethod pearson --clustMethod ward.D2 --dendroHeight 1 --percentDetection 90

# Load necessary packages (install them if it's not the case)
requiredPackages <- c('optparse',
                     'psych',
                     'ape')

for (p in requiredPackages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, repos = "https://cran.univ-paris1.fr/")
    
  }
  
  suppressMessages(invisible(library(p, character.only = TRUE)))
  
}

# Arguments
option_list <- list(
  make_option(
    "--file",
    default = NA,
    type = 'character',
    help = "Input file that contains expression values in log2(CPM+1). Columns = Observations (Samples, i.e. cells), Rows = Variables (Genes)."
  ),
  
  make_option(
    "--file_sep",
    default = "\t",
    type = 'character',
    help = "Column separator for the expression data file [default : '%default' ]"
  ),
  
  make_option(
    "--analysis_name",
    default = NA,
    type = 'character',
    help = "Name given to the current analysis"
  ),
  
  make_option(
    "--dataset",
    default = NA,
    type = 'character',
    help = "Details about the dataset under study that should appear in the output file names. Example = TestDataset_log2CPMvalues"
  ),
  
  make_option(
    "--corrMethod",
    default = "pearson",
    type = 'character',
    help = "Method to use for the correlation analysis between genes [default : '%default' ]"
  ),
  
  make_option(
   "--clustMethod",
    default = "ward.D2",
    type = 'character',
    help = "Method to use for the clustering analysis of genes [default : '%default' ]"
  ),
  
  make_option(
    "--dendroHeight",
    default = 1,
    type = 'double',
    help = "Height at which to cut the dendrogram [default : '%default' ]"
  ),
  
  make_option(
    "--percentDetection",
    default = 90,
    type = 'double',
    help = "Accepted percentage detection for HKG [default : '%default' ]"
  ),
  
  make_option(
    "--out",
    default = "./",
    type = 'character',
    help = "Path to set the working directory [default : '%default' ]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

# Define function to save csv files
save_csv <- function(dataframe_to_save, OutputNameFile) {
  write.csv(dataframe_to_save,
            file = paste(prefix,
                         OutputNameFile,
                         ".csv",
                         sep = "_"))
}

# Import normalized expression data
data.normCounts <- read.table(
  opt$file,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$file_sep,
  check.names = FALSE,
  row.names = 1
)

# Prefix for output files
prefix <- paste(Sys.Date(),  opt$dataset, opt$analysis_name, sep = "_")

# Setting Working Directory
dir.create(paste(opt$out, prefix, sep = "_"))

setwd(dir = paste(opt$out, prefix, sep = "_"))

# Compute descriptive statistics for each gene, including mean expression, standard deviation, coefficient of variation (CV), and detection frequency.
GeneMetadata <- data.frame(
  Min = apply(data.normCounts, 1, min), 
  Median = apply(data.normCounts, 1, median), 
  Mean = apply(data.normCounts, 1, mean), 
  SD = apply(data.normCounts, 1, sd), 
  Var = apply(data.normCounts, 1, var),
  Max = apply(data.normCounts, 1, max)
)  

GeneMetadata$Coeff_Var <- (GeneMetadata$SD/ GeneMetadata$Mean)*100 ## Coefficient of variation

GeneMetadata$NonDetection <- apply(data.normCounts, 1, function(x) sum(x ==0))  ## Number of cells in which the gene is not detected

GeneMetadata$Detection <- ncol(data.normCounts) - GeneMetadata$NonDetection  ## Number of cells in which the gene is detected

GeneMetadata$Percentage_NonDetection <- (GeneMetadata$NonDetection / ncol(data.normCounts)) *100  ## Percentage of cells in which the gene is not detected

GeneMetadata$Percentage_Detection <- 100 - GeneMetadata$Percentage_NonDetection  ## Percentage of cells in which the gene is detected

save_csv(GeneMetadata, "DescriptiveStatistics_Genes")

GeneMetadata$gene <- rownames(GeneMetadata)

# Compute the gene correlation matrix by calculating pairwise gene correlations across all cells.
corr <- cor(as.data.frame(t(data.normCounts)), method= opt$corrMethod)

# Convert the correlation matrix into a dissimilarity matrix
DissimilarityMatrix <- as.dist(1-corr)

# Perform hierarchical clustering of genes
## Cluster genes according to their correlation profiles

res.hc <- hclust(DissimilarityMatrix, method = opt$clustMethod)

pdf(file = paste0(prefix, "Dendrogram of Genes.pdf"))  ## Plot the obtained dendrogram

plot(res.hc, cex = 0.6, hang = -1)

dev.off()

# Define gene clusters

## Partition the hierarchical tree into gene clusters by cutting the dendrogram at the specified height.
### h is first set to the default value (=1).  Once the hierarchical tree is generated, we look at the hierarchical tree and decide at which height to cut it.
### The dendrogram cut height (--dendroHeight) directly influences the number and composition of gene clusters and may require optimization for different datasets.
CutTree <- cutree(res.hc, h=opt$dendroHeight)

NumberGenesPerCluster <- as.data.frame(table(CutTree))

save_csv(NumberGenesPerCluster, "Number of genes per cluster")

## Assign each gene to its corresponding cluster
GeneClusters <- data.frame(Cluster = CutTree)
GeneClusters$gene <- rownames(GeneClusters)

save_csv(GeneClusters, "Gene Clusters")

pdf(file = paste0(Sys.Date(), prefix, "Dendrogram of Genes_h=", opt$dendroHeight, ".pdf")) 

plot(res.hc, cex = 0.6, hang = -1)

rect.hclust(res.hc, h=opt$dendroHeight)

dev.off()

# Merge cluster assignments with the gene metadata
GeneMetadata <- merge(x= GeneClusters, y=GeneMetadata,by = "gene")

rownames(GeneMetadata) <- GeneMetadata$gene

save_csv(GeneMetadata, "GeneMetadata")

# Identify the least variable gene cluster
## Calculate the mean coefficient of variation per gene group
GeneGroups <- unique(GeneMetadata$Cluster) # a vector for unique group numbers

Variation <- as.data.frame(GeneMetadata$Coeff_Var)

rownames(Variation) <- rownames(GeneMetadata)

ClusterGenes <- lapply(GeneGroups, function(x) subset(Variation, GeneMetadata$Cluster==x))

MeanVariationPerGeneGroup <- as.data.frame(lapply(ClusterGenes, function(x) colMeans(x)))

colnames(MeanVariationPerGeneGroup) <- GeneGroups

rownames(MeanVariationPerGeneGroup) <- "AverageCoeffVariation"

MeanVariationPerGeneGroup <- as.data.frame(t(MeanVariationPerGeneGroup))

MeanVariationPerGeneGroup$Cluster <- rownames(MeanVariationPerGeneGroup)

MeanVariationPerGeneGroup  <- MeanVariationPerGeneGroup[order(MeanVariationPerGeneGroup$Cluster),] 

save_csv(MeanVariationPerGeneGroup, "MeanVariationPerGeneCluster")

## Select the cluster exhibiting the lowest average coefficient of variation
### The candidate housekeeping gene cluster is defined as the cluster exhibiting the lowest average coefficient of variation across all genes
cutoff <- MeanVariationPerGeneGroup[MeanVariationPerGeneGroup == min(MeanVariationPerGeneGroup$AverageCoeffVariation),2]

## Select candidate housekeeping genes
CandidateHKGgroup <- subset(GeneMetadata, GeneMetadata$Cluster==cutoff)

print(range(CandidateHKGgroup$Percentage_Detection))

print(range(CandidateHKGgroup$Coeff_Var))

save_csv(CandidateHKGgroup, paste0("CandidateHKGgroup", cutoff))

CandidateHKGgroup_percentZerosSelection <- subset(CandidateHKGgroup, CandidateHKGgroup$Percentage_Detection >= opt$percentDetection)

CandidateHKGgroup_percentZerosSelection <-  as.data.frame(rownames(CandidateHKGgroup_percentZerosSelection))

colnames(CandidateHKGgroup_percentZerosSelection) <- "gene"

save_csv(CandidateHKGgroup_percentZerosSelection, paste0("CorrGenes_LeastVar_DetectedIn", opt$percentDetection, "percentCells"))

save.image(file = paste0(prefix, ".RData"))

## The procedure identifies a set of candidate housekeeping genes that satisfy three criteria:
### (1) High detection frequency across single cells,
### (2) Low expression variability across the dataset, and
### (3) Membership in a highly correlated gene co-expression module.
