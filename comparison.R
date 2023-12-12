library(dplyr)
library(ggplot2)
library(VennDiagram)

edgeR_120 <- read.csv("edgeR_120_vs_control.csv")
edgeR_30 <- read.csv("edgeR_30_vs_control.csv")
deseq2_120 <- read.csv("contrast_120_vs_control.csv", row.names = 1)
deseq2_30 <- read.csv("contrast_30_vs_control.csv", row.names = 1)

# Filter for differentially expressed genes (DEGs) using adjusted p-value < 0.05
edgeR_120_DEGs <- subset(edgeR_120, PValue < 0.05)
edgeR_30_DEGs <- subset(edgeR_30, PValue < 0.05)
deseq2_120_DEGs <- subset(deseq2_120, padj < 0.05)
deseq2_30_DEGs <- subset(deseq2_30, padj < 0.05)

# Convert rownames to a column for DESeq2 data
deseq2_120_DEGs$Gene <- rownames(deseq2_120_DEGs)
deseq2_30_DEGs$Gene <- rownames(deseq2_30_DEGs)

# Identify common and unique genes
common_genes_120 <- intersect(edgeR_120_DEGs$Gene, deseq2_120_DEGs$Gene)
common_genes_30 <- intersect(edgeR_30_DEGs$Gene, deseq2_30_DEGs$Gene)

# Function to create and save a Venn diagram
save_venn_diagram <- function(gene_lists, file_name) {
  png(file_name, width = 800, height = 800) # Set the size of the image
  venn.plot <- venn.diagram(
    x = gene_lists,
    filename = NULL,
    output = TRUE
  )
  grid.draw(venn.plot)
  dev.off()
}

# Create and save Venn diagrams for the 120-minute and 30-minute comparisons
save_venn_diagram(list(edgeR = edgeR_120_DEGs$Gene, DESeq2 = deseq2_120_DEGs$Gene), "venn_120.png")
save_venn_diagram(list(edgeR = edgeR_30_DEGs$Gene, DESeq2 = deseq2_30_DEGs$Gene), "venn_30.png")


# Scatter plots for log2 fold change comparison
# Match genes in both datasets
matched_120 <- inner_join(edgeR_120_DEGs, deseq2_120_DEGs, by = "Gene")
matched_30 <- inner_join(edgeR_30_DEGs, deseq2_30_DEGs, by = "Gene")

# Scatter plot for 120 minutes
png("scatter_120.png")
ggplot(matched_120, aes(x = logFC, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(x = "edgeR Log2 Fold Change", y = "DESeq2 Log2 Fold Change", title = "120 minutes comparison")
dev.off()

# Scatter plot for 30 minutes
png("scatter_30.png")
ggplot(matched_30, aes(x = logFC, y = log2FoldChange)) +
  geom_point(alpha = 0.5) +
  theme_minimal() +
  labs(x = "edgeR Log2 Fold Change", y = "DESeq2 Log2 Fold Change", title = "30 minutes comparison")
dev.off()

# Correlation of log2 fold changes
cor_120 <- cor.test(matched_120$logFC, matched_120$log2FoldChange)
cor_30 <- cor.test(matched_30$logFC, matched_30$log2FoldChange)

print(paste("Correlation for 120 minutes: ", cor_120$estimate))
print(paste("Correlation for 30 minutes: ", cor_30$estimate))
