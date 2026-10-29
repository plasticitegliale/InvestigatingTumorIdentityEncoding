###########################################################
#   Identify genes specifically detected in a patient     #
###########################################################

# define arguments
Variables <- list(
  
  # Path to input file containing normalized expression values. Columns = Observations (Samples, i.e. cells), Rows = Variables (Genes).
  file= "NormalizedExpressionData_OutputStep3.tsv",
  
  # Column separator of input file
  sep = "\t", 
  
  # Name given to the current analysis
  analysis_name = "CommonGenesAllTumors",
  
  # Details about the dataset under study that should appear in the output file names
  dataset = "DatasetTest",
  
  # Path to input file containing cells’ metadata
  metadata = "CellMetadata.csv", 
  
  # Column separator for metadata file
  delimiter = ",", 

  # Percentage of patients in which a gene has to be detected to retain it
  percentDetection = 100,
  
  # Name of column that contains identity of tumor of origin
  tumor = "Sample.name", 
  
  # Path to set the working directory
  out = "."
)

# Load necessary packages (install them if it's not the case)
requiredPackages <- c('venn')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
  }
  suppressMessages(invisible(library(p, character.only = TRUE)))
}

# Import input files
  data.normCounts <- read.table(
    Variables$file,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = Variables$sep,
    check.names = FALSE,
    row.names = 1
  ) ## Columns = Observations (Samples), Rows = Variables (Genes)
  
  metadata <- read.csv(
    Variables$metadata,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = Variables$delimiter,
    check.names = FALSE,
    row.names = 1
  ) ## Rows = Observations (Samples), Columns = Variables

  metadata <- metadata[colnames(data.normCounts),]
  metadata$CellID <- rownames(metadata)

# Identifying the genes detected in each tumor
# A gene is considered detected if its expression is nonzero.
  detected <- data.normCounts != 0
  
## Tumor identities
  samples <- unique(as.character(metadata[[Variables$tumor]])) # a vector for unique group numbers

## Number of cells in each tumor
  cells_per_tumor <- table(
    factor(metadata[[Variables$tumor]], levels = samples)
  )

## Number of cells in which each gene is detected, per tumor
  Genes_detected_A <- sapply(samples, function(sample_id) {
    sample_cells <- rownames(metadata)[as.character(metadata[[Variables$tumor]]) == sample_id]
    rowSums(detected[, sample_cells, drop = FALSE], na.rm = FALSE)
  })      
  write.csv(Genes_detected_A, file=paste(Sys.Date(), "_Number of cells per tumor in which each gene is detected.csv")) 

# Define number of samples in which a gene has to be detected for it to be retained
  sample.thre <- (Variables$percentDetection * length(samples))/100

# Identify genes detected in every tumor.
## A gene is considered as detected in a tumor if detected in at least 1 cell from that tumor
  genes_detected_per_tumor <- Genes_detected_A > 0

  CommonDetectedGenes <- rownames(Genes_detected_A)[rowSums(genes_detected_per_tumor) >= sample.thre]

  write.csv(data.frame(gene = CommonDetectedGenes), file = file.path(Variables$out, paste0(Sys.Date(), "_", Variables$dataset, "_", Variables$analysis_name, "_CommonGenes.csv")), row.names = FALSE)

# Extract the expression matrix containing genes detected in all tumors.
  data.normCounts <- data.normCounts[rownames(data.normCounts) %in% CommonDetectedGenes$gene,]
  
  write.csv(data.normCounts, file = paste0(Sys.Date(), "_", Variables$dataset, "_ExpressionDataOf_", Variables$analysis_name, ".csv")) 
