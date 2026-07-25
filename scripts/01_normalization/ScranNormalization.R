#############################################################################
#              Normalize data using scran package per metacell              #
#############################################################################

# This Rscript allows to normalize data by a size factor calculated across cells or per metacell.
# The input file is the cells' expression data in count values.
# Output files include the scran-normalized expression matrix, which can be used for grouping analyses.

# Example of command (that generates output files) :
# Rscript ScranNormalization.R --file ExpressionData.csv --file_sep , --absolute_detection 3 --metacellNorm FALSE --out NormalizedExpressionData_OutputStep5_

# Load necessary packages (install them if it's not the case)
requiredPackages <- c('optparse', 'scran', 'scater')

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
    help = "Input file that contains expression values in counts. Columns = Observations (Samples, i.e. cells), Rows = Variables (Genes)."
  ),
  
  make_option(
    "--file_sep",
    default = "\t",
    type = 'character',
    help = "Column separator for the expression data file [default : '%default' ]"
  ),
  
  make_option(
    "--absolute_detection",
    default = 0,
    type = 'numeric',
    help = "Include genes with detected expression in at least this number of cells [default : '%default' ]"
  ),
  
  make_option(
    "--metacellNorm",
    default = FALSE,
    type = 'logical',
    help = "Whether to calculate the size factor per metacell (if TRUE) or across all cells (if FALSE) [default : '%default' ]"
  ),
  
  make_option(
    "--clustMethod",
    default = "hclust",
    type = 'character',
    help = "The clustering method to use for grouping cells into metacells, which can be either hclust or igraph [default : '%default' ]"
  ),
  
  make_option(
    "--out",
    default = "~",
    type = 'character',
    help = "Path to set the working directory [default : '%default' ]"
  )
)

opt <- parse_args(OptionParser(option_list = option_list), args = commandArgs(trailingOnly = TRUE))

# Import gene expression data file
data.counts <- read.table(
  opt$file,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$file_sep,
  check.names = FALSE,
  row.names = 1
)

# Remove lowly detected genes
## This filtering step reduces the influence of sparsely detected genes during normalization.
kept_genes <- rowSums(data.counts != 0) >= opt$absolute_detection
data.counts <- data.counts[kept_genes, ]

# Construct a SingleCellExperiment object to store all data
sce <- SingleCellExperiment(list(counts = as.matrix(data.counts)))

# Determine whether size factor should be calculated per metacell or per cell
if (opt$metacellNorm == TRUE) {
  ### Metacell-based normalization generally improves size-factor estimation in heterogeneous datasets by reducing the impact of biological variability on normalization.
  
  ## Cluster cells to create metacells (optional but recommended)
  clusters <- quickCluster(sce, min.size = 100, method = opt$clustMethod)
  
  ## Sum factors
  sce <- computeSumFactors(sce, cluster = clusters)
  
  ## To be added to the output file name
  outputFileName <- paste0("larger.logcounts_scran_", opt$clustMethod)
  
} else {
  ## Sum factors
  sce <- computeSumFactors(sce)
  
  ## To be added to the output file name
  outputFileName <- "logcounts_scran"
  
}

# Normalize counts using the estimated size factors and compute log-transformed normalized expression values.
sce <- scater::normalize(sce)

# Export the normalized expression matrix
logcounts <- data.frame(genes = rownames(sce),
             logcounts(sce),
             check.names = FALSE)

write.table(
  logcounts,
  paste0(opt$out, outputFileName, ".tsv"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)