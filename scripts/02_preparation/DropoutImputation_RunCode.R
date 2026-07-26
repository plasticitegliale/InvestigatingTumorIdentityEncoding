##############################################
#              Dropout imputation            #
##############################################

# This Rscript allows to impute dropout events using the scImpute approach.
# Input file includes the cells' expression data in log2(CPM/TPM +1) values.
# Output file include the imputed gene expression matrix.

# Rcode ran on 4.2.0

# define arguments
Variables <- list(
  # Path to input file containing normalized expression values. Columns = Observations (Samples, i.e. cells), Rows = Variables (Genes).
  file= "NormalizedExpressionData_OutputStep3.tsv",
  
  # Column separator of input file
  file_sep = "\t",
  
  # Name given to the current analysis
  analysis_name = "DropoutImputation",
  
  # Details about the dataset under study that should appear in the output file names
  dataset = "DatasetTest",
  
  # Path to set the working directory
  out = "."
)

# Load necessary packages (install them if it's not the case)
requiredPackages <- c('scImpute','kernlab', 'parallel', 'penalized', 'doParallel', 'foreach', 'iterators','rsvd')

for (p in requiredPackages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, repos = "https://cran.univ-paris1.fr/")
  }
  suppressMessages(invisible(library(p, character.only = TRUE)))
}

# Define function to save csv files
save_csv <- function(dataframe_to_save, OutputNameFile) {
  write.csv(dataframe_to_save,
            file = paste(prefix, OutputNameFile, ".csv", sep = "_"))
}

# import functions for scImpute
source("DropoutImputation_functions.R")

# Prefix for output files
prefix <- paste(Sys.Date(),  Variables$dataset, Variables$analysis_name, sep = "_")

# Import gene expression matrix
count_lnorm <- read.table(file = Variables$file, sep = Variables$file_sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)

# Remove rows with NA
row.has.na <- apply(count_lnorm, 1, function(x){any(is.na(x))})
sum(row.has.na)
count_lnorm <- count_lnorm[!row.has.na,]

# Replace all zeros by a pseudo-value (= log10(1.01))
count_lnorm[count_lnorm == 0] <- log10(1.01)

# Run imputation of dropout events
res_imp <- imputation_model8(count = count_lnorm, labeled = FALSE, 
                             point = log10(1.01), drop_thre = 0.75,
                             Kcluster = 1,
                             out_dir = Variables$out, ncores = 10)

# Retrieve results
expr_imputed <- res_imp$count_imp
outliers <- res_imp$outlier

# Replace all values equal to log10(1.01) by 0
## these values were considered as zero values and thus, were not imputed.
expr_imputed[expr_imputed == log10(1.01)] <- 0
save_csv(expr_imputed, "ImputedExpressionMatrix")
