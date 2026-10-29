##############################
#      Data Normalization    #
##############################

# Transform raw counts into log2(CPM +1)

# Example of command
# Rscript DataNormalization.R --data ExpressionData.csv --sep "," --out NormalizedExpressionData_OutputStep3.tsv

# Load necessary packages (install them if it's not the case)
requiredPackages = c('optparse')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE, quietly = T)) {
    install.packages(p)
  }
  suppressPackageStartupMessages(suppressMessages(library(p, character.only = TRUE)))
}

#Arguments
option_list = list(
  make_option(
    c("-d", "--data"),
    default = NA,
    type = 'character',
    help = "Input file that contains count values to transform : Rows = genes and Columns = Samples"
  ),
  make_option(
    c("-s", "--sep"),
    default = '\t',
    type = 'character',
    help = "File separator [default : '%default' ]"
  ),
  make_option(
    c("-c", "--colnames"),
    default = TRUE,
    type = 'logical',
    help = "Consider first line as header ? [default : '%default' ]"
  ),
  make_option(
    c("-o", "--out"),
    default = "res.tab",
    type = 'character',
    help = "Output name [default : '%default' ]"
  )
)


opt <- parse_args(OptionParser(option_list = option_list),
                 args = commandArgs(trailingOnly = TRUE))

if (opt$data == "" & !(opt$help)) {
  stop("At least one argument must be supplied (count data).\n",
       call. = FALSE)
}

# Open files
data.counts <- read.table(
  opt$data,
  h = opt$colnames,
  row.names = 1,
  sep = opt$sep,
  check.names = F
)

# Convert raw counts to CPM
cpm <- t(t(data.counts) / colSums(data.counts)) * 1000000

# Convert CPM to log2(CPM +1)
log2cpm1p <- log2(cpm +1)

# Save file
write.table(
  log2cpm1p,
  paste(opt$out),
  col.names = opt$colnames,
  row.names = T,
  quote = F,
  sep = "\t"
)
