###############################################
#              Define Metacells               #
###############################################

# This Rscript allows to define small groups of cells with a minimum of 3 cells and a maximum of 65, designated here as metacells.
# Input files include the expression data in log2(CPM+1) and the cells' metadata file.
# Output files include the clustering results as well as a .Rdata file.

## Rscript ran on 4.2.0

# Example of command (that generates output files) :
# Rscript DefiningMetacells.R --file NormalizedExpressionData_OutputStep3.tsv --metadata CellMetadata.csv --metadata_sep "," --analysis_name MetacellIdentification --dataset DatasetTest --corrMethod pearson --clustMethod ward.D2 --tumorColor black,green,orange,pink

# Load necessary packages (install them if it's not the case)
requiredPackages <- c('optparse',
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
    "--metadata",
    default = NA,
    type = 'character',
    help = "File that contains cells' metadata. Rows = Cells, Columns = Variables."
  ),
  
  make_option(
    "--metadata_sep",
    default = "\t",
    type = 'character',
    help = "Column separator for the metadata file [default : '%default' ]"
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
    help = "Method to use for the correlation analysis between cells [default : '%default' ]"
  ),
  
  make_option(
    "--clustMethod",
    default = "ward.D2",
    type = 'character',
    help = "Method to use for the clustering analysis of cells [default : '%default' ]"
  ),
  
  make_option(
    "--tumor",
    default = "Sample.name",
    type = 'character',
    help = "Column name corresponding to tumor identity"
  ),
  
  make_option(
    "--tumorColor",
    default = NA,
    type = 'character',
    help = "Color to assign to each tumor. List of colors comma separated [default : '%default' ]"
  ),
  
  make_option(
    "--out",
    default = "~",
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

# Define function to cut the hierarchical tree at specific heights h
CutDendro <- function(Input, dendroHeight){
  CutTree <- cutree(res.hc, h=dendroHeight)
  Clusters <- data.frame(CutTree)
  Clusters$cell <- rownames(Clusters)
  NumberCellsPerCluster <- data.frame(table(CutTree))
  save_csv(Clusters, paste0("Cell Clusters_h=",dendroHeight))
  save_csv(NumberCellsPerCluster, paste0("Cell Clusters(frequency)_h=", dendroHeight, "A"))
}

# Import expression data file
data.normCounts <- read.table(
  opt$file,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$file_sep,
  check.names = FALSE,
  row.names = 1
) ## Columns = Observations (Samples), Rows = Variables (Genes)

# Import metadata file
metadata <- read.table(
  opt$metadata,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$metadata_sep,
  check.names = FALSE,
  row.names = 1
) ## Rows = Observations (Samples), Columns = Variables

# Reorder the metadata so that cells appear in the same order as in the expression matrix.
## Critical: Cell identifiers MUST be identical and in the same order in both files.
metadata <- metadata[colnames(data.normCounts),]

# Prefix for output files
prefix <- paste(Sys.Date(),  opt$dataset, opt$analysis_name, sep = "_")

# Setting Working Directory
dir.create(paste(opt$out, prefix, sep = "_"))

setwd(dir = paste(opt$out, prefix, sep = "_"))

# Calculate pairwise cell correlations
## Compute the correlation matrix between all cells.
corr <- cor(data.normCounts, method= opt$corrMethod)

## Convert the correlation matrix into a dissimilarity matrix.
DissimilarityMatrix <- as.dist(1-corr)

# Perform hierarchical clustering
## Cluster cells according to their transcriptomic similarity.
res.hc <- hclust(DissimilarityMatrix, method = opt$clustMethod)  

colorPalette  <- unlist(strsplit(opt$tumorColor, split = ","))

pdf(file = paste0(prefix, "Dendrogram of Cells.pdf"))  ## Plot the obtained dendrogram

plot(res.hc, cex = 0.2, hang = -1)

plot(as.phylo(res.hc), type = "fan", tip.color = colorPalette[as.factor(metadata[, opt$tumor])], label.offset = 0.0, cex = 0.2)

legend("topleft", legend=levels(as.factor(metadata[, opt$tumor])), cex= 0.75, col=colorPalette, pch=16)

dev.off()

# Define metacells
## Cut the hierarchical tree at specific heights based on dendrogram (to have groups of cells with min 3 cells and max 65 cells)
dendroHeight <- 1.1
CutTree <- cutree(res.hc, h = dendroHeight)
Clusters <- data.frame(Cell = names(CutTree), Metacell = CutTree)
NumberCellsPerCluster <- as.data.frame(table(CutTree))

## Repeat this procedure using multiple dendrogram heights to identify the partition that best satisfies the desired metacell size.
dendroHeight <- 1
CutTree <- cutree(res.hc, h = dendroHeight)
Clusters <- data.frame(Cell = names(CutTree), Metacell = CutTree)
NumberCellsPerCluster <- as.data.frame(table(CutTree))

dendroHeight <- 0.90
CutTree <- cutree(res.hc, h = dendroHeight)
Clusters <- data.frame(Cell = names(CutTree), Metacell = CutTree)
NumberCellsPerCluster <- as.data.frame(table(CutTree))

dendroHeight <- 0.98
CutTree <- cutree(res.hc, h = dendroHeight)
Clusters <- data.frame(Cell = names(CutTree), Metacell = CutTree)
NumberCellsPerCluster <- as.data.frame(table(CutTree))

# Refine metacell boundaries
## Inspect the cluster sizes obtained at the different dendrogram heights.
### The final metacell partition should contain groups comprising 3–65 transcriptionally similar cells, providing a balance between robustness and biological resolution.
### Groups containing more than 65 cells should be subdivided using a lower dendrogram height, whereas groups containing fewer than 3 cells should be merged with neighboring groups identified at a higher dendrogram height.
### Critical: Because dendrogram topology varies among datasets, the final metacell definition requires manual inspection to ensure that each metacell contains between 3 and 65 cells.

# Save analysis
write.csv(Clusters, "MetacellGroups.csv")
save.image(file = paste0(prefix, ".RData"))

# Critical notes
## Pearson correlation is used by default to quantify transcriptomic similarity between cells, although alternative correlation methods may be specified.
## The optimal dendrogram cut height is dataset-dependent and should be determined by evaluating the resulting metacell sizes.
## Metacells are intended to improve downstream analyses such as normalization and visualization and should not be interpreted as final biological cell populations.
