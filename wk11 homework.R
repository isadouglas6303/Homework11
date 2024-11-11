#nicePCA

library(ggplot2)

nicePCA <- function(pcaData, groupVec) {
  df <- data.frame(PC1 = pcaData$x[, 1],
                   PC2 = pcaData$x[, 2],
                   Group = as.character(groupVec))
  
  variance1 <- round(pcaData$pve[1] * 100, 2)
  variance2 <- round(pcaData$pve[2] * 100, 2)
  
  plot <- ggplot(df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point() +
    labs(x = paste("PCA 1 (", variance1, "% variance explained)", sep = ""),
         y = paste("PCA 2 (", variance2, "% variance explained)", sep = ""),
         title = "PCA Plot") +
    theme_minimal() +
    coord_fixed() 
  
  return(plot)
}


#niceVolcano

library(ggplot2)
library(ggrepel)

niceVolcano <- function(resData, pvalCutoff = 0.05, log2FCCutoff = 1.5) {
  # Check that 'geneName' exists
  if (!"geneName" %in% colnames(resData)) {
    resData$geneName <- rownames(resData)
  }
  
  # Remove rows with missing values
  resData <- resData[!is.na(resData$padj) & !is.na(resData$log2FoldChange), ]
  
  # Add columns for significance and cutoffs
  resData$Significant <- with(resData, padj < pvalCutoff)
  resData$ExceedsCutoff <- with(resData, abs(log2FoldChange) > log2FCCutoff)
  resData$PlotType <- with(resData, ifelse(Significant & ExceedsCutoff, "Significant & Exceeds Cutoff",
                                           ifelse(Significant, "Significant", "Not Significant")))
  
  # Volcano plot
  plot <- ggplot(resData, aes(x = log2FoldChange, y = -log10(padj), shape = PlotType, color = ExceedsCutoff)) +
    geom_point(size = 2) +
    scale_color_manual(values = c("grey", "turquoise")) +
    scale_shape_manual(values = c("Significant & Exceeds Cutoff" = 17, "Significant" = 16, "Not Significant" = 1)) +
    labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", 
         title = "Volcano Plot") +
    theme_minimal() +
    geom_hline(yintercept = -log10(pvalCutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-log2FCCutoff, log2FCCutoff), linetype = "dashed", color = "blue") +
    geom_text_repel(data = subset(resData, Significant & ExceedsCutoff),
                    aes(label = geneName), max.overlaps = 100)
  
  return(plot)
}


#resSummary

resSummary <- function(resData) {
  # Check that columns exist
  if (!"padj" %in% colnames(resData)) {
    stop("The input data frame must have a 'padj' column.")
  }
  if (!"log2FoldChange" %in% colnames(resData)) {
    stop("The input data frame must have a 'log2FoldChange' column.")
  }
  
  # Use 'gene_name' column for output
  if ("gene_name" %in% colnames(resData)) {
    resData$geneSymbol <- as.character(resData$gene_name)
  } else {
    resData$geneSymbol <- rownames(resData)
  }
  
  # Summary stats
  totalGenes <- nrow(resData)
  nonMissingPvals <- sum(!is.na(resData$padj))
  downregulated <- sum(resData$padj < 0.05 & resData$log2FoldChange < 0, na.rm = TRUE)
  upregulated <- sum(resData$padj < 0.05 & resData$log2FoldChange > 0, na.rm = TRUE)
  
  # Top 10 genes
  topGenes <- head(resData[order(resData$padj), "geneSymbol"], 10)
  
  # Print output
  message("How many genes were included in the DESeq2 analysis after filtering?")
  print(totalGenes)
  
  message("How many genes included in the DESeq2 analysis had non-missing p-adjusted values?")
  print(nonMissingPvals)
  
  message("Note: using 0.05 as the p-adjusted cutoff for significance for the following calculations.")
  
  message("How many genes were significantly downregulated in group2 relative to group1?")
  print(downregulated)
  
  message("How many genes were significantly upregulated in group2 relative to group1?")
  print(upregulated)
  
  message("Names of the top 10 most significantly differentially expressed genes from smallest to largest p-adjusted:")
  print(topGenes)
}

