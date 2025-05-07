
# RNA-Seq Differential Expression Analysis: Tumor vs Normal (GSE81089)
# Author: Ifthekar Hussain, Sri Vasista Talagampala
# Date: April 26, 2025

# Load necessary libraries
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# Load count data
counts <- read.table("GSE81089_readcounts_featurecounts1.tsv", header = TRUE, sep = "\t", row.names = 1)

# Load metadata
meta_lines <- readLines("GSE81089_series_matrix (2).txt")
char_line <- grep("!Sample_characteristics_ch1", meta_lines, value = TRUE)[1]
group <- ifelse(grepl("T\"?$", unlist(strsplit(char_line, "\t"))[-1]), "Tumor", "Normal")
sample_info <- data.frame(sample = colnames(counts), condition = factor(group))

# Clean the count data and prepare DESeq2 object
counts_clean <- counts[complete.cases(counts), ]
dds <- DESeqDataSetFromMatrix(countData = counts_clean, colData = sample_info, design = ~condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)
res <- results(dds)

# Volcano Plot
EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 1)

# Heatmap of top 30 DE genes
vsd <- vst(dds, blind = FALSE)
topGenes <- head(order(res$padj), 30)
mat <- assay(vsd)[topGenes, ]
annotation <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, annotation_col = annotation, scale = "row")

# PCA Plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot", x = paste0("PC1: ", percentVar[1], "%"), y = paste0("PC2: ", percentVar[2], "%"))

# Save significant genes
sig_genes <- na.omit(res[res$padj < 0.05, ])
write.csv(sig_genes, "significant_DEGs_GSE81089.csv")
