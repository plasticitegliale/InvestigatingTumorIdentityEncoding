#######################################################################################################
#              Normalize data by the expression of housekeeping genes (HKG) per metacell              #
#######################################################################################################

# This Rscript allows to normalize data by the expression of housekeeping genes (HKG). Because normalization factors are derived from groups of transcriptionally similar cells rather than individual cells, this strategy reduces technical variability while minimizing the influence of stochastic single-cell fluctuations.
# Input files include the cells' expression data in log2(CPM+1) values, and their metadata.
# Output files include the normalized expression matrices in log2(CPM+1), which can be used for grouping analyses.

## Rscript ran on 4.2.0

# Example of command (that generates output files) :
# Rscript DataNormalizationByHKG.R --file NormalizedExpressionData_OutputStep3.tsv --dataset DatasetTest --analysis_name Test_NormalizeByHKGPerMetacell --NormalizationGroup Metacell --HKG 2026-07-25_DatasetTest_HKGidentification_CandidateHKGgroup177_.csv --HKG_sep , --metadata MetacellGroups.csv --metadata_sep ,

# Load necessary packages (install them if it's not the case)
requiredPackages <- c('optparse',
                      'psych')

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
    help = "File that contains cells' metadata with metacell labels. Rows = Cells, Columns = Variables. File generated from 'DefiningMetacells.R' code." 
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
    "--NormalizationGroup",
    default = NA,
    type = 'character',
    help = "Group by which the data should be normalized"
  ),
  
  make_option(
    "--HKG",
    default = NA,
    type = 'character',
    help = "File containing the list of housekeeping genes. All genes in 1 column with gene as column names. One gene per row. Can be obtained by running the IdentifyingHKG.R code"
  ),
  
  make_option(
    "--HKG_sep",
    default = "\t",
    type = 'character',
    help = "Column separator for HKG file [default : '%default' ]"
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

# Define function to normalize gene expression values
NormByHKG_apply <- function(x) {
  apply(x, 2, function(y) (y / metadata[i, "NormalizingFactor"]))
}

# Import gene expression data file
data.normCounts <- read.table(
  opt$file,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$file_sep,
  check.names = FALSE,
  row.names = 1
) 

# Import metadata file
## Obtained after running 'DefiningMetacells.R'
metadata <- read.csv(
  opt$metadata,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$metadata_sep,
  check.names = FALSE,
  row.names = 1
)

# Critical: Cell identifiers must be identical and ordered consistently between the expression matrix and the metadata.
metadata <- metadata[colnames(data.normCounts), ]

# Import housekeeping gene (HKG) file
## Obtained after running 'IdentifyingHKG.R'
Genes <- read.csv(
  opt$HKG,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = opt$HKG_sep,
  check.names = FALSE,
  row.names = 1
) 

colnames(Genes) <- "gene"

# Prefix for output files
prefix <- paste(Sys.Date(),  opt$dataset, opt$analysis_name, sep = "_")

# Setting Working Directory
dir.create(paste(opt$out, prefix, sep = "_"))

setwd(dir = paste(opt$out, prefix, sep = "_"))

# Extract housekeeping gene expression
HKG_ExpressionData <- subset(data.normCounts, rownames(data.normCounts) %in% Genes$gene)

save_csv(HKG_ExpressionData, "HKG_ExpressionMatrix_log2(CPM+1)")

# Calculate the average housekeeping gene expression for each metacell

## Group cells according to their assigned metacells
metacells <- unique(metadata[,opt$NormalizationGroup]) # a vector for unique group numbers

## Calculate the average expression of each housekeeping gene within every metacell
CellsPerMetacell <- lapply(metacells, function(x)
  subset(t(HKG_ExpressionData), metadata[,opt$NormalizationGroup] == x))

MeanHKG_perMetacell <-  data.frame(lapply(CellsPerMetacell, function(x)
  colMeans(x)))

colnames(MeanHKG_perMetacell) <- metacells

# Compute metacell normalization factors
## Calculate the geometric mean of the average housekeeping gene expression within each metacell.
## The resulting geometric mean serves as the normalization factor for all cells belonging to the corresponding metacell.
NormalizingFactorPerMetacell <- data.frame(geometric.mean(MeanHKG_perMetacell))

colnames(NormalizingFactorPerMetacell) <- "NormalizingFactor"

NormalizingFactorPerMetacell$metacell <- rownames(NormalizingFactorPerMetacell)

save_csv(NormalizingFactorPerMetacell, "HKGGeomMeanPerMetacell")

# Are the HKG detected in all metacells?
## Critical: Housekeeping genes should be consistently detected across metacells. Metacells containing numerous undetected housekeeping genes may produce unreliable normalization factors.
MeanHKG_perMetacell$N_zeros <- apply(MeanHKG_perMetacell, 1, function(x)
  sum(x == 0))

save_csv(MeanHKG_perMetacell, "HKGMetadata")

# Assign a normalizing factor to each cell depending on which metacell it belongs to
metadata$NormalizingFactor <- metadata[,opt$NormalizationGroup]

for (i in 1:nlevels(as.factor(metadata[,opt$NormalizationGroup]))) {
  metadata$NormalizingFactor[metadata$NormalizingFactor == i] <- NormalizingFactorPerMetacell[i, 1]
}

metadata <- metadata[colnames(data.normCounts), ]

save_csv(metadata, "MetadataWithNormalizationFactor")

# Normalize each gene expression by the calculated normalisation factor
data.normCounts_inv <- as.data.frame(t(data.normCounts))

HKGnormalized_data.normCounts <- data.frame(matrix(, nrow = nrow(data.normCounts_inv), ncol = ncol(data.normCounts_inv)))

for (cell in 1:nrow(data.normCounts_inv)) {
  HKGnormalized_data.normCounts[cell, ] <- NormByHKG_apply(data.normCounts_inv[cell, ])
}

colnames(HKGnormalized_data.normCounts) <- colnames(data.normCounts_inv)

rownames(HKGnormalized_data.normCounts) <- rownames(data.normCounts_inv)

HKGnormalized_data.normCounts <- as.data.frame(t(HKGnormalized_data.normCounts))

save_csv(HKGnormalized_data.normCounts, "HKGnormalizedExpressionMatrix")

# save.image(file = paste0(prefix, ".RData"))