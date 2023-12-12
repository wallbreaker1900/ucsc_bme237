library(edgeR)
library(ggplot2)

count_files <- c(
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890279_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890278_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890277_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890276_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890275_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890274_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890273_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890272_HISAT2_sorted.bam_counts.txt",
  "/hb/home/tyu41/project/htseq_counts_output/SRR13890271_HISAT2_sorted.bam_counts.txt"
)

# Conditions for each sample, correctly ordered based on your description
conditions <- factor(c("control", "control", "control",
                       "treatment_30min", "treatment_30min", "treatment_30min",
                       "treatment_120min", "treatment_120min", "treatment_120min"))

# Read counts function
read_counts <- function(file) {
  read.table(file, header = FALSE, row.names = 1)
}

# Read and merge count data
count_data_list <- lapply(count_files, read_counts)

# Combine the count data with unique column names
count_matrix <- do.call(cbind, count_data_list)
colnames(count_matrix) <- paste("Sample", seq_along(count_files), sep = "_")

# edgeR analysis
dge <- DGEList(counts = count_matrix, group = conditions)
dge <- dge[filterByExpr(dge), ]
dge <- calcNormFactors(dge)
design <- model.matrix(~ conditions)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)

# Corrected glmQLFTest calls using the actual coefficient names from the design matrix
qlf_120_vs_control <- glmQLFTest(fit, coef = "conditionstreatment_120min")
qlf_30_vs_control <- glmQLFTest(fit, coef = "conditionstreatment_30min")

# Convert results to data frames and check
results_120_vs_control <- as.data.frame(qlf_120_vs_control$table)
results_120_vs_control$Gene <- rownames(dge$counts)

results_30_vs_control <- as.data.frame(qlf_30_vs_control$table)
results_30_vs_control$Gene <- rownames(dge$counts)

# Save results to CSV
write.csv(results_120_vs_control, file = "/hb/home/tyu41/project/edgeR_120_vs_control.csv", row.names = FALSE)
write.csv(results_30_vs_control, file = "/hb/home/tyu41/project/edgeR_30_vs_control.csv", row.names = FALSE)

create_volcano_plot_edgeR <- function(result_df, title, file_name) {
  # Apply Benjamini-Hochberg correction for the PValue column
  result_df$padj <- p.adjust(result_df$PValue, method = "BH")
  result_df$color <- NA
  
  # Assign color only to significant points based on padj and logFC
  result_df$color[result_df$padj < 0.05 & abs(result_df$logFC) > 1] <- 'green'  # Highly upregulated
  result_df$color[result_df$padj < 0.05 & abs(result_df$logFC) > 0.6 & abs(result_df$logFC) <= 1] <- 'red'  # Moderately upregulated
  result_df$color[result_df$padj < 0.05 & abs(result_df$logFC) <= 0.6] <- 'blue'  # Significant but low magnitude

  # Filter out non-significant genes by removing rows with NA color
  result_df <- na.omit(result_df)

  # Create the plot for only significant points
  p <- ggplot(result_df, aes(x = logFC, y = -log10(PValue), color = color)) +
    geom_point(size = 1.5) +
    scale_color_identity() +
    theme_bw() +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    theme(legend.position = "none")  # If you do not wish to display a legend
    xlim(c(-3, 3))
  ggsave(file_name, plot = p, width = 6, height = 4, dpi = 300, path = "/hb/home/tyu41/project/")
}

create_volcano_plot_edgeR(results_120_vs_control, "Volcano Plot (120 min vs Control)", "volcano_120_vs_control_edgeR.png")
create_volcano_plot_edgeR(results_30_vs_control, "Volcano Plot (30 min vs Control)", "volcano_30_vs_control_edgeR.png")

padj_threshold <- 0.05
log2FC_up_threshold <- 0.6
log2FC_down_threshold <- -0.6

# Function to sort, filter, and save significantly regulated genes for edgeR results
save_sorted_differentially_expressed_genes_edgeR <- function(results_df, up_file_name, down_file_name) {
  results_df$padj <- p.adjust(results_df$PValue, method = "BH")
  results_df$abs_logFC <- abs(results_df$logFC)

  upregulated <- subset(results_df, padj < padj_threshold & logFC > log2FC_up_threshold)
  downregulated <- subset(results_df, padj < padj_threshold & logFC < log2FC_down_threshold)

  upregulated_sorted <- upregulated[order(-upregulated$abs_logFC),]
  downregulated_sorted <- downregulated[order(-downregulated$abs_logFC),]

  write.csv(upregulated_sorted, file = up_file_name, row.names = FALSE)
  write.csv(downregulated_sorted, file = down_file_name, row.names = FALSE)
}

# Apply the function for both conditions
save_sorted_differentially_expressed_genes_edgeR(results_120_vs_control, "/hb/home/tyu41/project/upregulated_120_vs_control_edgeR.csv", "/hb/home/tyu41/project/downregulated_120_vs_control_edgeR.csv")
save_sorted_differentially_expressed_genes_edgeR(results_30_vs_control, "/hb/home/tyu41/project/upregulated_30_vs_control_edgeR.csv", "/hb/home/tyu41/project/downregulated_30_vs_control_edgeR.csv")
