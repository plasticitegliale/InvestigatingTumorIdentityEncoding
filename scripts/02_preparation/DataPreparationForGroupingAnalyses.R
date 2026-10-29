#########################################################
#       Preparing data for cell grouping analysis       #
#########################################################

# Run prior to cell grouping analysis, this Rscript allows to:
#   A. exclude or retain some specific genes,
#   B. standardize data by tumor/patient,
#   C. binarize data,
#   D. randomly select different numbers of genes, or
#   E. identify most expressed genes
# Input files include the expression data in log2(CPM+1), the list of genes of interest and/or the cells' metadata file.
# Output files include the new expression matrices, which can be used as input for grouping analyses.

# Example of command (that generates output files) :
### Rscript DataPreparationForGroupingAnalyses.R --file NormalizedExpressionData_OutputStep3.tsv --file_sep "\t" --dataset DatasetTest --out . --analysis_name BinarizedData --binarize TRUE

# Load necessary packages (install them if not installed)
requiredPackages <- c('optparse', 'Hmisc')

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
    help = "Input file that contains normalized expression values. Columns = Observations (Cells), Rows = Variables (Genes)."
  ),
  
  make_option(
    "--file_sep",
    default = "/t",
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
    help = "Details about the dataset under study that should appear in the output file names. Example = TestDataset_log2CPMvalues."
  ),
  
  make_option(
    "--gene",
    default = NA,
    type = 'character',
    help = "File that contains the genes of interest (csv format). 1 column with column name = gene. Required if some genes should be excluded, or only those genes should be retained."
  ),
  
  make_option(
    "--gene_sep",
    default = "/t",
    type = 'character',
    help = "Column separator for the gene file [default : '%default' ]"
  ),
  
  make_option(
    "--exclude",
    default = FALSE,
    type = 'logical',
    help = "Genes from the gene file are excluded if TRUE. Else, only these genes are retained. [default : '%default' ]"
  ),
  
  make_option(
    "--groups",
    default = NA,
    type = 'character',
    help = "Data standardization will be done for each of these groups separately"
  ),
  
  make_option(
    "--metadata",
    default = NA,
    type = 'character',
    help = "File that contains cells' metadata. Required only if data is to be standardized by tumor/patient. Rows = Observations (Cells), Columns = Features on each cell (e.g. individual of origin, cell type, ...)."
  ),
  
  make_option(
    "--metadata_sep",
    default = "/t",
    type = 'character',
    help = "Column separator for the metadata file [default : '%default' ]"
  ),
  
  make_option(
    "--binarize",
    default = FALSE,
    type = 'logical',
    help = "Whether or not to binarize expression data [default : '%default' ]"
  ),
  
  make_option(
    "--NumberOfGenes",
    default = NA,
    type = 'integer',
    help = "The number of genes to randomly select"
  ),
  
  make_option(
    "--Replicate",
    default = NA,
    type = 'integer',
    help = "The number of replicates for which the random gene selection should be done"
  ),
  
  make_option(
    "--ExpressionCutoff",
    default = NA,
    type = 'integer',
    help = "The cutoff to use to identify most expressed genes"
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

#  Setting Working Directory
setwd(dir = opt$out)

# Define function to save csv files
save_csv <- function(dataframe_to_save, OutputNameFile, prefix_output = prefix) {
    write.csv(dataframe_to_save, file = paste(prefix_output, OutputNameFile, ".csv", sep = "_"))
  }

# Prefix for output files
prefix <- paste(Sys.Date(),  opt$dataset, opt$analysis_name, sep = "_")

# Import expression data file
data.normCounts <- read.table(
    opt$file,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = opt$file_sep,
    check.names = FALSE,
    row.names = 1
  )

# Retain specific genes only or exclude them
if(!is.na(opt$gene)) {
  ## Import the symbols of the genes to be excluded or retained
  Genes <- read.csv(
    opt$gene,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = opt$gene_sep,
    check.names = FALSE,
    row.names = 1
  ) 
  
  genes <- as.character(Genes$gene)
  
  if (opt$exclude == TRUE) {
    ## Exclude the expression data of the genes in the Genes dataframe
    Filtered_Dataset <-
      data.normCounts[which(rownames(data.normCounts) %nin% genes),]
    
    save_csv(Filtered_Dataset, "ExpressionMatrix")
    
  } else {
    ## Retain only the expression data of the genes in the Genes dataframe
    Filtered_Dataset <- data.normCounts[!(rownames(data.normCounts) %nin% genes),]
    
    save_csv(Filtered_Dataset, "ExpressionMatrix")
  }
}

# Standardize expression data by tumor/patient
if (!is.na(opt$groups)) {
  ## Import the cells' metadata file
  metadata <- read.csv(
    opt$metadata,
    header = TRUE,
    stringsAsFactors = FALSE,
    sep = opt$metadata_sep,
    check.names = FALSE,
    row.names = 1
  )
  
  ## To make sure that the order of the cells in the expression data file and in the metadata file are the same
  metadata <- metadata[colnames(data.normCounts), ]
  
  ## Get the individual IDs
  IndividualNames <- as.list(levels(as.factor(metadata[, opt$groups])))
  
  ## Standardize data per tumor/patient
  StandardizedData <- lapply(IndividualNames, function(individual_index) {
    Individual <- as.data.frame(subset(t(data.normCounts), metadata[, opt$groups] == individual_index))
    
    data_standardized_for_one_tumor <- as.data.frame(t(scale(Individual)))
    
    return(data_standardized_for_one_tumor)
  })
  
  ## Combine the individual standardized data into a single dataframe
  StandardizedData_dataframe <- do.call("cbind", StandardizedData)
  
  ## If standard deviation = 0, the standardized value = NA. Therefore, all NA are replaced by 0.
  StandardizedData_dataframe <- replace(StandardizedData_dataframe, is.na(StandardizedData_dataframe), 0)
  
  save_csv(StandardizedData_dataframe, "StandardizedExpressionMatrix")
  
}

# Binarize expression data
if (opt$binarize == TRUE) {
  ## Binarize into an ON-OFF system.
  ## For each cell, detected genes were considered as ON and assigned a value of 1 regardless of their relative expression levels.
  ## Undetected genes were considered as OFF and retained their original value of 0.
  data.normCounts[data.normCounts > 0] <- 1
  
  save_csv(data.normCounts, "BinarizedExpressionMatrix")
  
}

# Extract the expression values, across all cells, of different numbers of randomly-selected genes  
if (!is.na(opt$NumberOfGenes)) {
  lapply(1:opt$Replicate, function(n_replicates) {
    set.seed(123)
    
    RandomGenes_Dataset <- data.normCounts[sample(nrow(data.normCounts), opt$NumberOfGenes), ]
    
    save_csv(RandomGenes_Dataset, paste0("Replicate", n_replicates, "_ExpressionMatrix"))
    }
  )
}

# Select most expressed genes
if (!is.na(opt$ExpressionCutoff)) {
  ## prefix for output files
  prefix_expGenes_analysis <- paste(
    Sys.Date(),
    opt$dataset,
    paste0("TopMostExpGenes_ExpressionCutoff=", opt$ExpressionCutoff),
    sep = "_"
  ) 
  
  ## Calculate the CPM values from the log2(CPM+1) values
  data.CPM <- (2 ^ data.normCounts) - 1
  
  ## Calculate the mean CPM values for each gene
  data.CPM$Mean <- rowMeans(data.CPM)
  
  ## Calculate the log2(meanCPM) values for each gene
  data.CPM$logMeanCPM <- log2((data.CPM$Mean) + 1)
  
  ## Select genes whose mean expression is above the cutoff, and extract their expression values across cells
  MostExpressed_data.CPM <- data.CPM[data.CPM$logMeanCPM >= opt$ExpressionCutoff, ]
  
  ## Exclude the columns "Mean" and "logMeanCPM" from the expression matrix
  MostExpressed_data.CPM <-
    subset(MostExpressed_data.CPM, select = !(
      colnames(MostExpressed_data.CPM) %in% c("Mean", "logMeanCPM")
    ))
  
  ## Calculate the log2(CPM+1) values from the CPM values
  MostExpressed_data.logCPM <- log2((MostExpressed_data.CPM) + 1)
  
  save_csv(MostExpressed_data.logCPM,
           paste0("ExpresssionMatrix_log2CPM"),
           prefix_output = prefix_expGenes_analysis)
  
}
