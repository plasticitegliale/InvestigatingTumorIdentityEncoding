######################################################################################################################
#              Grouping analysis based on Hierarchical Clustering on Principal Components (HCPC) approach            #
######################################################################################################################

# This Rscript allows to perform unsupervised grouping analysis using the HCPC approach.

# R code ran on R 4.2.0

# Input files include:
## the cells' normalized expression data (for e.g. in log2(CPM+1) values)
## the cells' metadata

# Output files include the clustering results, a .RData file as well as some plots (PCA/tSNE/UMAP, heatmap).

# Define arguments
Variables <- list(
  # Path to input file containing normalized expression values. Columns = Observations (Samples, i.e. cells), Rows = Variables (Genes).
  file= "NormalizedExpressionData_OutputStep3.tsv",
  
  # Column separator of input file
  file_sep = "\t",
  
  # Path to input file containing cells' metadata. Rows = Cells, Columns = Variables.
  metadata = "CellMetadata.csv",
  
  # Column separator for the metadata file
  metadata_sep = ",",
  
  # Name given to the current analysis. This name will be added to the output file names.
  analysis_name = "NoDataManipulation",
  
  # Details about the dataset under study that should appear in the output file names.
  dataset = "DatasetTest",
  
  # Number of principal components to consider for grouping analyses
  numberPCsToUseForClustering = 10,
  
  # Number of groups in which to group cells. If set to NA, optimal cluster number as determined by the Silhouette index will be used.
  nClusters = NA, 
  
  # Column name corresponding to tumor identity
  tumor =  "Sample.name",

  # Color to assign to each tumor. List of comma-separated colors.
  tumorColor = "black,green,orange,pink", 
  
  # Whether or  not to run tSNE for data visualization
  tsne = FALSE,
  
  # Whether or  not to use UMAP for data visualization
  UMAP = TRUE,
  
  # Path to set the working directory
  out = "."
)

# Load necessary packages (install them if it's not the case)
requiredPackages <- c(
  'FactoMineR',
  'factoextra',
  'NbClust',
  'RColorBrewer',
  'tsne',
  'ClusterR',
  'ggplot2',
  'circlize',
  'ComplexHeatmap',
  'umap',
  'psych',
  'RSpectra'
)

for (p in requiredPackages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, repos = "https://cran.univ-paris1.fr/")
  }
  suppressMessages(invisible(library(p, character.only = TRUE)))
}

# Import required functions
source("CellGroupingAnalysis_HCPC_functions.R")

# Import expression data file
data.normCounts <- read.table(
  Variables$file,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = Variables$file_sep,
  check.names = FALSE
)

# Import metadata file
metadata <- read.csv(
  Variables$metadata,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = Variables$metadata_sep,
  check.names = FALSE,
  row.names = 1)

metadata <- metadata[colnames(data.normCounts),]

# Prefix for output files
prefix <- paste(Sys.Date(),  Variables$dataset, Variables$analysis_name, sep = "_")

# Set working directory
setwd(Variables$out)

# Run PCA
res.pca <- PCA(
  t(data.normCounts),
  ncp = Variables$numberPCsToUseForClustering,
  graph = FALSE,
  scale.unit = FALSE
)

Cells_coord_onPCA <- as.data.frame(res.pca$ind$coord)

# Extracting the variance data from the PCA
Vartab <- get_eig(res.pca)

# Assess the clustering tendency of the dataset using Hopkins statistic (H)
## Function to calculate Hopkins statistic in the package: H = sum(minp)/(sum(minp) + sum(minq))
## If H > 0.5, the dataset is not uniformly distributed, i.e., it contains meaningful clusters.
## If H < 0.5, the dataset is uniformly distributed, i.e., there is no meaningful clusters.

clustend <- get_clust_tendency(Cells_coord_onPCA, seed = 123, n = (nrow(Cells_coord_onPCA) - 1))
clustend$hopkins_stat

# Identify cell groups using HCPC
if (clustend$hopkins_stat > 0.5) {
  if(is.na(Variables$nClusters) == TRUE){
    
    ## Determining the optimal number of clusters using Silhouette method
    OptimalClusterNumber_Silhouette <- NbClust(data = Cells_coord_onPCA, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 30, method = "ward.D2", index = "silhouette")
    Ncluster <- as.numeric(OptimalClusterNumber_Silhouette$Best.nc[1])
    
  }else{
    ## Manually set number of clusters based on dendrogram
    Ncluster <- Variables$nClusters 
  }
  
  # Perform Hierarchical Clustering On Principal Components Followed By Kmeans Clustering (HCPC).
  pdf(file = paste0(prefix, "_Dendrogram.pdf"))
  res.hcpc <- HCPC(res.pca, nb.clust = Ncluster, graph = TRUE)   
  dev.off()
  
  # Saving Grouping Results
  save_csv(res.hcpc$desc.var$quanti.var, "QuanVarDescribingClusters")
  save_multiple_files(as.list(res.hcpc$desc.var$quanti), "QuanVarDescribingCluster")
  save_csv(res.hcpc$desc.axes$quanti.var,"AxesDescribingClusters")
  save_multiple_files(as.list(res.hcpc$desc.axes$quanti), "AxesDescribingCluster")
  
  # Update metadata file
  metadata$Cluster <- res.hcpc$data.clust[, (nrow(data.normCounts) + 1)]
  save_csv(metadata, "MetadataWithClusterIDs")

 # Save table with number of cells per tumor or per cluster
  save_csv(table(metadata[, Variables$tumor]), "tumor_sizes")
  save_csv(table(metadata$Cluster), " cluster_sizes")

  # NMI score calculation
  sink(paste0(Sys.Date(), "_External-validation-ContributionOfEachTumorToEachCluster.txt"))
  external_validation(as.numeric(factor(metadata[, Variables$tumor])),
                      as.numeric(metadata$Cluster), summary_stats = TRUE)
  sink()
  
} else {
  stop(
    paste("Hopkins statistic (H) =", clustend$hopkins_stat,   "(Since H < 0.5, the dataset is uniformly distributed, i.e., there is no meaningful clusters.")
  )
}

# Save outputs from PCA (scree plot and cell groups on PCA space)
if (clustend$hopkins_stat > 0.5) {
  
  pdf(file = paste0(prefix, "_PCA_scree plot.pdf"))
  fviz_eig(res.pca) ## Scree Plot
  dev.off()
  
  pdf(file = paste0(prefix, "_PCAplots.pdf"))
  plot_function(SampleCoord = Cells_coord_onPCA, AssignBasedOn = "tumor", PlotType = "PCA") ## Color by tumor of origin
  plot_function(SampleCoord = Cells_coord_onPCA, AssignBasedOn = "cluster", PlotType = "PCA") ## Color by cluster
  dev.off()
}

# Highlight cell groups on tSNE space
## tSNE visualization of HCPC clustering
if (clustend$hopkins_stat > 0.5) {
  if (Variables$tsne == TRUE) {
    set.seed(1) # Set random seed for reproducibility
    
    # define function to monitor embedding along tSNE optimization
epoch_callback <- function(embedding) {
  plot(
    embedding, pch = 16, col = "black", main = "t-SNE optimization progress"
   )
}    

  ## Calculate pairwise Euclidean distances between cells
    DistanceMatrix <- dist(Cells_coord_onPCA, method = "euclidean")
    
  ## Perform t-SNE dimensionality reduction
    Cells_coord_onTSNE <- tsne(DistanceMatrix, initial_config = NULL, k = 10, initial_dims = 30, perplexity = 50, max_iter = 1000, min_cost = 0, epoch_callback = epoch_callback, whiten = TRUE, epoch = 100)
    
    Cells_coord_onTSNE <- as.data.frame(Cells_coord_onTSNE)
    rownames(Cells_coord_onTSNE) <- rownames(Cells_coord_onPCA)
    save_csv(Cells_coord_onTSNE, "SampleCoordinates_OnTSNEplot")
    
    pdf(file = paste(prefix,  "tSNEPlots.pdf", sep = "_"))
    plot_function(SampleCoord = Cells_coord_onTSNE, AssignBasedOn = "tumor", PlotType = "tSNE")
    plot_function(SampleCoord = Cells_coord_onTSNE, AssignBasedOn = "cluster", PlotType = "tSNE")
    dev.off()
  }
}

# Highlight cell groups on UMAP space
if (clustend$hopkins_stat > 0.5) {
  if (Variables$UMAP == TRUE) {
    set.seed(123)
    Data.umap <- umap(d = Cells_coord_onPCA, init = "spectral")
    layout <- Data.umap
    if (class(Data.umap) == "umap") {
      layout <- Data.umap$layout
    }
    save_csv(layout, "UMAP_SampleCoordinates")
    pdf(file = paste(prefix,  "UMAPplots.pdf", sep = "_"))
    plot_function(SampleCoord = layout, AssignBasedOn = "tumor", PlotType = "UMAP")
    plot_function(SampleCoord = layout, AssignBasedOn = "cluster", PlotType = "UMAP")
    dev.off()
  }
}

# Save the main results, input parameters and R package versions used
## Save the input parameters used
sink(paste0(prefix, "_InputParameters.txt"))
cat(
  '\nUsing :\n- Data file (--file) :', Variables$file,
  '\n- Metadata file containing tumor of origin or other cell features (--metadata) :', Variables$metadata,
  '\n- OutputFolder (--out) :', Variables$out,
  '\n\nOther parameters :\n- Type of analysis (--analysis_name) :', Variables$analysis_name,
  '\n- Name of dataset (--dataset) :', Variables$dataset,
  '\n- Number of principal components considered for grouping analyses (--numberPCsToUseForClustering) :', Variables$numberPCsToUseForClustering,
  '\n- Column name corresponding to tumor identity (--tumor) :', Variables$tumor,
  '\n- Color assigned to each tumor (--tumorColor) :', Variables$tumorColor,
  '\n- Should tSNE be run (TRUE = Yes and FALSE  = No) (--tsne) :', Variables$tsne,
  '\n- Should UMAP be run (TRUE = Yes and FALSE  = No) (--UMAP) :', Variables$UMAP,
  '\n'
)
sink()

## Save the versions of R and R packages used for reproducibility
sink(paste0(prefix, "_RpackagesVersions.txt"))
sessionInfo()
sink()

## Save the output reporting standards
sink(paste0(prefix, "_DatasetInformation.txt"))
cat('\nTotal number of cells in dataset :', nrow(metadata),
  '\n- Number of genes analyzed :', nrow(data.normCounts), 
  '\n- Number of tumors :', length(unique(metadata[, Variables$tumor])),
  '\n- Number of clusters :', length(unique(metadata$Cluster)),
  '\n'
)
sink()





