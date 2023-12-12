library(DESeq2)
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

# Function to read counts and return as a data frame
read_counts <- function(file) {
  read.table(file, header = TRUE, row.names = 1)
}

# Read and merge count data
count_data_list <- lapply(count_files, read_counts)
count_matrix <- do.call(cbind, count_data_list)

# Ensure that rownames (gene names) are properly set
rownames(count_matrix) <- rownames(count_data_list[[1]])

# Create a sample information data frame
sampleInfo <- data.frame(
  sampleName = sub(".txt", "", basename(count_files)),
  condition = c("control", "control", "control",
                "treatment_30min", "treatment_30min", "treatment_30min",
                "treatment_120min", "treatment_120min", "treatment_120min")
)

# Create DESeqDataSet from the matrix
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sampleInfo,
                              design = ~ condition)

# Perform DESeq2 analysis
dds <- DESeq(dds)

# Define the two comparisons of interest
contrast_120_vs_control <- results(dds, contrast = c("condition", "treatment_120min", "control"))
contrast_30_vs_control <- results(dds, contrast = c("condition", "treatment_30min", "control"))

# Save the results to CSV files
write.csv(contrast_120_vs_control, "contrast_120_vs_control.csv", row.names = TRUE)
write.csv(contrast_30_vs_control, "contrast_30_vs_control.csv", row.names = TRUE)

create_volcano_plot <- function(contrast_data, title) {
  # Convert DESeqResults to a dataframe
  contrast_df <- as.data.frame(contrast_data)
  
  # Apply Benjamini-Hochberg correction
  contrast_df$padj <- p.adjust(contrast_df$pvalue, method = "BH")

  # Assign colors based on fold change and adjusted p-value
  contrast_df$color <- ifelse(contrast_df$padj < 0.05 & abs(contrast_df$log2FoldChange) <= 0.6, "Not DE",
                              ifelse(contrast_df$padj < 0.05 & contrast_df$log2FoldChange > 1, "Highly Upregulated",
                                     ifelse(contrast_df$padj < 0.05 & contrast_df$log2FoldChange > 0.6, "Moderately Upregulated", NA)))

  # Filter out NA (non-significant) colors
  contrast_df <- contrast_df[!is.na(contrast_df$color), ]

  ggplot(contrast_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = color), size = 1.5) +
    scale_color_manual(values = c("Not DE" = "blue", 
                                  "Moderately Upregulated" = "red", 
                                  "Highly Upregulated" = "green")) +
    theme_bw() + 
    labs(
      x = "Log2 Fold Change",
      y = "-log10(p-value)",
      title = title
    ) +
    guides(color = guide_legend(title = NULL)) + 
    xlim(c(-3, 3))
}

volcano_plot_120_vs_control <- create_volcano_plot(contrast_120_vs_control, "Volcano Plot (120 min vs Control)")
ggsave("volcano_120_vs_control.png", plot = volcano_plot_120_vs_control, width = 6, height = 4)

volcano_plot_30_vs_control <- create_volcano_plot(contrast_30_vs_control, "Volcano Plot (30 min vs Control)")
ggsave("volcano_30_vs_control.png", plot = volcano_plot_30_vs_control, width = 6, height = 4)

padj_threshold <- 0.05
log2FC_up_threshold <- 0.6
log2FC_down_threshold <- -0.6

# Function to sort, filter, and save significantly regulated genes
save_sorted_differentially_expressed_genes <- function(contrast, up_file_name, down_file_name) {
  contrast$abs_log2FC <- abs(contrast$log2FoldChange)
  upregulated <- subset(contrast, padj < padj_threshold & log2FoldChange > log2FC_up_threshold)
  downregulated <- subset(contrast, padj < padj_threshold & log2FoldChange < log2FC_down_threshold)

  upregulated_sorted <- upregulated[order(-upregulated$abs_log2FC),]
  downregulated_sorted <- downregulated[order(-downregulated$abs_log2FC),]

  write.csv(upregulated_sorted, file = up_file_name, row.names = TRUE)
  write.csv(downregulated_sorted, file = down_file_name, row.names = TRUE)
}

# Apply the function for both contrasts
save_sorted_differentially_expressed_genes(contrast_120_vs_control, "upregulated_120_vs_control.csv", "downregulated_120_vs_control.csv")
save_sorted_differentially_expressed_genes(contrast_30_vs_control, "upregulated_30_vs_control.csv", "downregulated_30_vs_control.csv")
