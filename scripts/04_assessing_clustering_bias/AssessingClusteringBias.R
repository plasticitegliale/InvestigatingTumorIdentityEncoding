##########################################################################################################
#              Assessing whether cell clusters are enriched for cells from a specific feature            #
##########################################################################################################

# This Rscript allows to evaluate the contribution of different cell features to the identified cell clusters.
# R code ran on R 4.2.0

# Input files include clustering results and cells' annotations. Clustering results can be obtained from the CellGroupingAnalysis_HCPC_RunCode.R

# Output files include the normalized mutual information (NMI) score, other clustering metrics as well as an alluvial plot as visualization.

# Define variables
Variables <- list(
  
  ## path to clustering results and cell features. Table in csv format with cells as rows and each feature as column
  cluster.path = "2026-07-25_DatasetTest_NoDataManipulation_MetadataWithClusterIDs_.csv", 
  
  ## name of column containing cell groups (from clustering analysis)
  cell.group.col = "Cluster", 
  
  ## name of column containing the cell feature to assess. Cell feature should be a qualitative label
  feature.compare = "Location", # Sample.name, Location
  
  ## name to add to output file
  out.prefix = "Region", # Tumor, Region
  
  # Color to assign to each group. List of comma-separated colors.
  Color = "black,green,orange,pink"
  
)

# import cells’ annotations
metadata <- read.csv(Variables$cluster.path, header=TRUE, row.names=1)

# Visualize clustering bias
## prepare data
library(dplyr)
data <- table(metadata[, Variables$feature.compare], metadata[, Variables$cell.group.col]) %>% as.data.frame()
names(data) <- c(Variables$feature.compare, Variables$cell.group.col, "CellNumber")

## plot
library(ggplot2)
library(ggalluvial)
colorPalette  <- unlist(strsplit(Variables$Color, split = ","))
pdf(paste0(Sys.Date(), "_AlluvialPlot_ContributionOfEach", Variables$out.prefix, "ToEachCluster.pdf"), width = 3, height = 3)
ggplot(data,
       aes(axis1 = data[,Variables$cell.group.col],
           axis2 = data[,Variables$feature.compare],
           y = CellNumber)) +
  geom_alluvium(aes(fill = data[,Variables$feature.compare]), alpha = 0.8, width = 0.5) +
  geom_stratum(width = 0.5) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void() +
  scale_fill_manual(values = colorPalette) +
  theme(legend.position = "bottom") +
  labs(fill = Variables$feature.compare)
dev.off()

# Quantify clustering bias
## Compute the normalized mutual information (NMI) between cell feature and the assigned cell clusters
library(ClusterR)
sink(paste0(Sys.Date(), "_External-validation-ContributionOfEach", Variables$out.prefix, "ToEachCluster.txt"))

external_validation(
  as.numeric(factor(metadata[, Variables$feature.compare])),
  as.numeric(factor(metadata[, Variables$cell.group.col])),
  summary_stats = TRUE
)

sink()

# Define function to evaluate cluster–tumor correspondence using NMI
label_permutation <- function(cluster.ids, tumor.ids, n_perm = 10000, seed = 123) {
  
  # Calculate observed NMI
  observed_nmi <- external_validation(as.numeric(factor(tumor.ids)), as.numeric(factor(cluster.ids)), method = "nmi")
  
  # Generate the null distribution by permuting tumor labels
  set.seed(seed)
  null_nmi <- replicate(n_perm, {
    permuted_labels <- sample(tumor.ids, replace = FALSE)
    external_validation(as.numeric(factor(permuted_labels)), as.numeric(factor(cluster.ids)), method = "nmi")
  })
  
  # Calculate permutation p-value
  p_value <- (1 + sum(null_nmi >= observed_nmi)) / (n_perm + 1)
  
  # Summarize results as a table
  results <- data.frame(
    Observed_NMI = observed_nmi,
    Mean_Null_NMI = mean(null_nmi),
    SD_Null_NMI = sd(null_nmi),
    NMI_Enrichment = observed_nmi - mean(null_nmi),
    Permutation_P_value = p_value,
    Number_of_Permutations = n_perm
  )
  return(list(summary = results, null_distribution = null_nmi))
}

# Perform label permutation
results <- label_permutation(cluster.ids = metadata[, Variables$cell.group.col], tumor.ids = metadata[, Variables$feature.compare], n_perm = 10000, seed = 123)
write.csv(results$summary, paste0(Sys.Date(), "_ObservedNMIrelativeToNullNMIdistribution.csv"))
NMI.null <- as.data.frame(results$null_distribution)
names(NMI.null) <- "NMIval"

# Visualize data
library(ggplot2)
pdf(paste0(Sys.Date(), "_Histogram_ObservedNMIrelativeToNullNMIdistribution.pdf"))
ggplot(NMI.null, aes(x=NMIval)) + 
  geom_histogram(bins = 40) +
  theme_classic() +
  xlab("NMI under permuted tumor labels") + ylab("Frequency") + ggtitle("NMI permutation test") +
  labs(caption = paste0("ObsNMI=", round(results$summary$Observed_NMI, 4), " ; pval=", format.pval(NMI.res$Permutation_P_value, 4))) +
  scale_x_log10() +
  geom_vline(xintercept = results$summary$Observed_NMI, linetype = "dashed", color = "red")
dev.off()

# Define function for repeated balanced subsampling
cell_subsampling <- function(cell_metadata, cluster_col, patient_col, label_col, n_per_patient = NULL, n_repeats = 500, seed = 123) {
  
  # Extract relevant data from cell annotation table
  cell_metadata_subset <- cell_metadata[, c(cluster_col, patient_col, label_col)]
  names(cell_metadata_subset) <- c("cluster", "patient", "label")
  
  # Get number of cells from each patient
  patient_counts <- table(cell_metadata_subset$patient)
  if (length(patient_counts) < 2) {
    stop("At least two patients are required.")
  }
  
  # Determine the number of cells to sample per patient
  if (is.null(n_per_patient)) {
    n_per_patient <- min(patient_counts)
  }
  
  if (n_per_patient < 1 || n_per_patient != as.integer(n_per_patient)) {
    stop("n_per_patient must be a positive integer.")
  }
  
  if (n_per_patient > min(patient_counts)) {
    stop("n_per_patient exceeds the number of cells available for at least one patient. Reduce n_per_patient or exclude patients using a prespecified criterion.")
  }
  
  # Calculate original observed NMI
  original_nmi <- external_validation(as.numeric(factor(cell_metadata_subset$label)), as.numeric(factor(cell_metadata_subset$cluster)), method = "nmi")
  
  # Set random seed for reproducibility
  set.seed(seed)
  
  # Create list to store NMI for every subsampling
  subsampling_results <- vector("list", n_repeats)
  
  for (i in seq_len(n_repeats)) {
    # Sample the same number of cells from each patient
    sampled_cells <- do.call(
      rbind,
      lapply(names(patient_counts), function(patient_id) {
        
        patient_cells <- cell_metadata_subset[cell_metadata_subset$patient == patient_id,, drop = FALSE]
        selected_rows <- sample(seq_len(nrow(patient_cells)), size = n_per_patient, replace = FALSE)
        patient_cells[selected_rows, , drop = FALSE]
      })
    )
    
    # Calculate NMI for the sampled cells
    nmi_value <- external_validation(as.numeric(factor(sampled_cells$label)), as.numeric(factor(sampled_cells$cluster)), method = "nmi")
    
    subsampling_results[[i]] <- data.frame(
      repetition = i,
      NMI = nmi_value,
      n_cells = nrow(sampled_cells),
      n_patients = length(unique(sampled_cells$patient))
    )
  }
  subsampling_results <- do.call(rbind, subsampling_results)
  
  # Summarize variability across subsamplings
  summary_results <- data.frame(
    Original_NMI = original_nmi,
    Sampled_NMI_Median = median(subsampling_results$NMI),
    Sampled_NMI_Mean = mean(subsampling_results$NMI),
    Sampled_NMI_SD = sd(subsampling_results$NMI),
    Sampled_NMI_IQR = IQR(subsampling_results$NMI),
    Percentile_2.5 = unname(quantile(subsampling_results$NMI, 0.025)),
    Percentile_97.5 = unname(quantile(subsampling_results$NMI, 0.975)),
    N_patients = length(patient_counts),
    Cells_per_patient = n_per_patient,
    Total_sampled_cells = length(patient_counts) * n_per_patient,
    N_repeats = n_repeats
  )
  
  # Return results
  return(list(summary = summary_results, repetitions = subsampling_results, original_NMI = original_nmi))
}

# Perform balanced subsampling test
results <- cell_subsampling(
  cell_metadata = metadata,
  cluster_col = "Cluster",
  patient_col = "Sample.name",
  label_col = "Sample.name",
  n_per_patient = NULL,  # if NULL, uses the smallest patient cell count
  n_repeats = 500,
  seed = 123
)

write.csv(results$summary, paste0(Sys.Date(), "_ObservedNMIrelativeToNMIfromsubsampling_Summary.csv"))
write.csv(results$repetitions, paste0(Sys.Date(), "_ObservedNMIrelativeToNMIfromsubsampling_Repetitions.csv"))

# Visualize data
library(ggplot2)
pdf(paste0(Sys.Date(), "_Histogram_ObservedNMIrelativeToNMIfromsubsampling.pdf"))
ggplot(results$repetitions, aes(x=NMI)) + 
  geom_histogram(bins = 40) +
  theme_classic() +
  xlab("NMI value") +
  ylab("Frequency") +
  ggtitle("NMI robustness under balanced subsampling") +
  labs(caption = paste0("OriginalNMI=", round(results$summary$Original_NMI, 4))) +
  scale_x_log10() +
  geom_vline(xintercept = results$summary$Original_NMI, linetype = "dashed", color = "red")

results$summary$grp <- "SubsamplingNMI"
ggplot(results$summary, aes(x=grp, y = Sampled_NMI_Mean)) +
  geom_col(alpha = 0.8, color = "black") +
  geom_point(aes(y = Original_NMI), shape = 22, size = 5, fill = "gold", color = "black") +
  geom_errorbar(aes(ymin = Sampled_NMI_Mean - Sampled_NMI_SD, ymax = Sampled_NMI_Mean + Sampled_NMI_SD), width = 0.2, color = "black", linewidth = 0.8) +
  labs(title = "NMI robustness under balanced subsampling", x = "", y = "NMI distribution") +
  theme_classic()
dev.off()


