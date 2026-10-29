### Identifying highly variable genes using scran ###
#####################################################

# Note: In comparison with other methods (BASiCS, Brennecke, scLVM, scVEGs and Seurat), scran has a good performance in clustering, stable number of detected HVGs, good independency from the mean and good running time. (https://doi.org/10.1093/bib/bby011)

# Define arguments
Variables <- list(
  
  # Path to input file containing normalized expression values. Columns = Observations (Samples, i.e. cells), Rows = Variables (Genes).
  file= "NormalizedExpressionData_OutputStep3.tsv",
  
  # Column separator of input file
  file_sep = "\t",
  
  # Number of highly variable genes to retain or exclude
  nHVG = 500,
  
  # Whether to keep or exclude the top highly variable genes. Can be “keep” or “exclude”
  filter = "keep",
  
  # Name of the output file 
  output_name = "Top500HVG"
)

# Load necessary packages
library(scran)
library(Hmisc)

# Import expression data file
data.normCounts <- read.table(file = Variables$file, sep=Variables$file_sep, header = TRUE, row.names = 1)

# Identify top highly variable genes
detectVarGenes <- function(mat, ngenes=500, FDR=0.00001){
  var.out = scran::modelGeneVar(mat)

  hvg.out = var.out[which(var.out$FDR <= FDR),]
  hvg.out = hvg.out[order(hvg.out$bio, decreasing=TRUE),]
  genes=rownames(hvg.out)[1:ngenes]
  return(genes)
}

# Detect highly variable genes
hvg <- detectVarGenes(mat=as.matrix(data.normCounts), ngenes=Variables$nHVG, FDR = 1e-05)

# Check the number of detected HVGs
nHVG_available <- length(hvg)
if (nHVG_available == 0) {
  stop("No significant highly variable genes were detected at FDR <= 1e-5. Check the input data or reconsider the FDR threshold.")
}else{
if (nHVG_available < Variables$nHVG) {
  warning(paste0("Requested ", Variables$nHVG, " highly variable genes, but only ", nHVG_available, " were returned by detectVarGenes(). Using all available genes.")))
     }
}

# Filter expression matrix
if(Variables$filter == "keep"){
  Filtered_Dataset <- data.normCounts[hvg,]
} else {
  Filtered_Dataset <- data.normCounts[which(rownames(data.normCounts) %nin% hvg),]
}
write.csv(Filtered_Dataset, file= paste0(Sys.Date(), Variables$output_name, "ExpressionMatrix.csv"))
