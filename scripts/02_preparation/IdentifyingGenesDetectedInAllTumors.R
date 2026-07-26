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
samples <- unique(metadata[, Variables$tumor]) # a vector for unique group numbers
SampleCells <- lapply(samples, function(x) subset(t(data.normCounts), metadata[, Variables$tumor] ==x))

Genes_detected <- lapply(SampleCells, function(x) (colSums(x != 0)))
Genes_detected_A <- as.data.frame(Genes_detected)
colnames(Genes_detected_A) <- samples
write.csv(Genes_detected_A, file=paste(Sys.Date(), "_Number of cells per tumor in which each gene is detected.csv"))

# Compare list of genes detected in each tumor
CompareGenes <- list()
for(i in 1:ncol(Genes_detected_A)){
  CompareGenes[i] <- list(as.character(rownames(Genes_detected_A[Genes_detected_A[,i] != 0,])))
}

if(length(CompareGenes)>7){
  nCompareLists <- length(CompareGenes) / 7
  if(nCompareLists >1 & nCompareLists <= 2){
    CompareGenes.venn1 <- venn(CompareGenes[1:7], simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection1 <- attr(CompareGenes.venn1, "intersections")
    
    CompareGenes.venn2 <- venn(CompareGenes[8:length(CompareGenes)], simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection2 <- attr(CompareGenes.venn2, "intersections")
    
    CompareLists <- list()
    CompareLists["A"] <- Intersection1[length(Intersection1)]
    CompareLists["B"] <- Intersection2[length(Intersection2)]
    CompareGenes.venn <- venn(CompareLists, simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection <- attr(CompareGenes.venn, "intersections")
  }
  
  if(nCompareLists >2 & nCompareLists <= 3){
    CompareGenes.venn1 <- venn(CompareGenes[1:7], simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection1 <- attr(CompareGenes.venn1, "intersections")
    
    CompareGenes.venn2 <- venn(CompareGenes[8:14], simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection2 <- attr(CompareGenes.venn2, "intersections")
    
    CompareGenes.venn3 <- venn(CompareGenes[15:length(CompareGenes)], simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection3 <- attr(CompareGenes.venn3, "intersections")
    
    CompareLists <- list()
    CompareLists["A"] <- Intersection1[length(Intersection1)]
    CompareLists["B"] <- Intersection2[length(Intersection2)]
    CompareLists["C"] <- Intersection3[length(Intersection3)]
    CompareGenes.venn <- venn(CompareLists, simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
    Intersection <- attr(CompareGenes.venn, "intersections")
  }
} else{
  CompareGenes.venn <- venn(CompareGenes, simplify = FALSE, intersections = TRUE, show.plot=TRUE, zcolor = "style")
  Intersection <- attr(CompareGenes.venn, "intersections")
}

# Get list of genes detected in all tumors
CommonDetectedGenes <- as.data.frame(Intersection[length(Intersection)])
colnames(CommonDetectedGenes) <- "gene"

# Extract expression values of genes detected in all tumors
data.normCounts <- data.normCounts[rownames(data.normCounts) %in% CommonDetectedGenes$gene,]
write.csv(data.normCounts, file = paste0(Sys.Date(), "_", Variables$dataset, "_ExpressionDataOf_", Variables$analysis_name, ".csv"))
