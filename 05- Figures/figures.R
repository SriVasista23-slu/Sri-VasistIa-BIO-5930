# Load required packages
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Assumes 'dds' and 'res' already exist in the environment

# 1. PCA Plot
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  labs(
    title = "PCA Plot",
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )
ggsave("figures/pca_plot.png", p)

# 2. Volcano Plot
png("figures/volcano_plot.png", width = 800, height = 600)
EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1,
    title = "Volcano Plot: Tumor vs Normal")
dev.off()

# 3. MA Plot
ma_data <- as.data.frame(res)
ma_data$significance <- ifelse(is.na(ma_data$padj), "NA",
                          ifelse(ma_data$padj < 0.05, "Significant", "Not Significant"))
p2 <- ggplot(ma_data, aes(x = baseMean, y = log2FoldChange, color = significance)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  labs(title = "MA Plot")
ggsave("figures/ma_plot.png", p2)

# 4. Heatmap
topGenes <- head(order(res$padj), 30)
mat <- assay(vsd)[topGenes, ]
annotation <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
pheatmap(mat, annotation_col = annotation, scale = "row",
         main = "Heatmap of Top 30 Genes",
         filename = "figures/heatmap_top_genes.png")

# 5. Top gene boxplot
top_gene <- rownames(res)[which.min(res$padj)]
plot_data <- plotCounts(dds, gene = top_gene, intgroup = "condition", returnData = TRUE)
p3 <- ggplot(plot_data, aes(x = condition, y = count)) +
  geom_jitter(width = 0.2, color = "blue") +
  geom_boxplot(alpha = 0.3) +
  labs(title = paste("Expression of", top_gene))
ggsave("figures/top_gene_count_plot.png", p3)

# 6. Log2FC density
p4 <- ggplot(res, aes(x = log2FoldChange)) +
  geom_density(fill = "skyblue", alpha = 0.6) +
  labs(title = "Density of log2 Fold Changes")
ggsave("figures/log2fc_density.png", p4)

# 7. Adjusted p-value histogram
p5 <- ggplot(res, aes(x = padj)) +
  geom_histogram(binwidth = 0.01, fill = "gray") +
  labs(title = "Histogram of Adjusted P-values")
ggsave("figures/padj_histogram.png", p5)

# 8. Top 10 genes barplot
res_clean <- na.omit(res)
res_sorted <- res_clean[order(res_clean$log2FoldChange), ]
top_genes <- rbind(head(res_sorted, 10), tail(res_sorted, 10))
p6 <- ggplot(top_genes, aes(x = reorder(rownames(top_genes), log2FoldChange), 
                            y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_col() + coord_flip() +
  labs(title = "Top Up/Down Regulated Genes", x = "Gene", y = "log2 Fold Change")
ggsave("figures/top_genes_barplot.png", p6)
