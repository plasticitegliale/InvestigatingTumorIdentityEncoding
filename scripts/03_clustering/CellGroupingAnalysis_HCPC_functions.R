####################################################################################################################################
#              Functions for Grouping analysis based on Hierarchical Clustering on Principal Components (HCPC) approach            #
####################################################################################################################################

# This Rscript contains functions needed to perform unsupervised grouping analysis using the HCPC approach.

# Descriptive Statistics Function
descriptive_stats <- function(InputData, OutputFileName) {
  SummaryData <- data.frame(
    mean = rowMeans(InputData),
    SD = apply(InputData, 1, sd),
    Variance = apply(InputData, 1, var),
    Percentage_Detection = apply(InputData, 1, function(x, y = InputData) {
      (sum(x != 0) / ncol(y)) * 100
    })
  )
  
  write.csv(SummaryData, file = paste0(Sys.Date(), "_DescriptiveStats_", OutputFileName, ".csv"))
  
  return(SummaryData)
}

# Define function to save csv files
save_csv <- function(dataframe_to_save, OutputNameFile) {
  write.csv(dataframe_to_save,
            file = paste(prefix,
                         OutputNameFile,
                         ".csv",
                         sep = "_"))
}

# Define function to save multiple csv files at once
save_multiple_files <- function(list_to_save, OutputNameFile) {
  i <- 1
  
  while (i <= length(list_to_save)) {
    write.csv(list_to_save[i],
              file = paste(prefix,
                           OutputNameFile,
                           i,
                           ".csv",
                           sep = "_"))
    
    i <- i + 1
  }
}

# Define function for 2D PCA or tSNE plots
plot_function <- function(SampleCoord, AssignBasedOn, PlotType) {
  ## Assign a color to each cell on the plot
  if (AssignBasedOn == "tumor") {
    ### Get the tumor names
    SampleLegend  <- levels(as.factor(metadata[, Variables$tumor]))
    
    ### Define the color palette to be used to color cells
    colorPalette  <- unlist(strsplit(Variables$tumorColor, split = ","))
    
    ### Assign a color to each cell based on tumor of origin
    SampleColor  <- colorPalette[as.factor(metadata[, Variables$tumor])]
    
  } else {
    ### Get the cluster names
    SampleLegend  <- levels(as.factor(metadata$Cluster))
    
    ### Define the color palette to be used to color cells
    qual_col_pals  <- brewer.pal.info[brewer.pal.info$category  == 'qual', ]
    
    colorPalette  <- unique(unlist(mapply(
      brewer.pal,
      qual_col_pals$maxcolors,
      rownames(qual_col_pals)
    )))
    
    ### Assign a color to each cell based on clustering results
    SampleColor <- vector() ## initiate variable
    for (cell_index in 1:length(metadata$Cluster)) {
      SampleColor[cell_index]  <- colorPalette[which(SampleLegend == metadata$Cluster[cell_index])]
    }
  }
  
  ## Define all axes to be plotted
  AxisCombinationsList <-
    lapply(1:(ncol(SampleCoord) - 1), function(xaxis) {
      PossibleXaxis <- xaxis
      
      yaxisPlot <- (xaxis + 1):ncol(SampleCoord)
      
      Combinations <- expand.grid(PossibleXaxis, yaxisPlot)
    })
  
  AxisCombinations_dataframe <- do.call("rbind", AxisCombinationsList)
  
  ## Plot the XY graph
  lapply(1:nrow(AxisCombinations_dataframe), function(CombinationNumber) {
    xaxisPlot <- AxisCombinations_dataframe[CombinationNumber, 1]
    
    yaxisPlot <- AxisCombinations_dataframe[CombinationNumber, 2]
    
    ### Define the label for the x-axis and y-axis
    if (PlotType == "PCA") {
      xlabel <- paste("PC", xaxisPlot, "(", format(round(Vartab[xaxisPlot, 2], 2)), " %)")
      
      ylabel <- paste("PC", yaxisPlot, "(", format(round(Vartab[yaxisPlot, 2], 2)), " %)")
      
    } else {
      xlabel <- paste("Component", xaxisPlot)
      
      ylabel <- paste("Component", yaxisPlot)
    }
    
    ### XY plot
    plot(
      SampleCoord[, c(xaxisPlot, yaxisPlot)],
      type = "p",
      pch = 16,
      col = SampleColor,
      cex = 1.5,
      xlab = xlabel,
      ylab = ylabel,
      axes = TRUE
    )
    
    #### Add legend on the XY plot
    legend(
      "topleft",
      legend = SampleLegend,
      col = colorPalette,
      pch = 16,
      cex = 1.5
    )
    
    ### empty return because otherwise legend coordinates are returned
    return()
    
    ### XY plot without legend
    plot(
      SampleCoord[, c(xaxisPlot, yaxisPlot)],
      type = "p",
      pch = 16,
      col = SampleColor,
      cex = 1.5,
      xlab = xlabel,
      ylab = ylabel,
      axes = TRUE
    )
    
    
  })
}