# Figures

This folder contains all the visual outputs generated from the RNA-Seq differential gene expression analysis of lung tumor vs. normal tissues using dataset GSE81089.

Each figure helps interpret different aspects of the gene expression data.

##  Contents

| File Name                      | Description |
|-------------------------------|-------------|
| `pca_plot.png`                | PCA plot showing sample clustering between tumor and normal groups |
| `volcano_plot.png`           | Volcano plot showing significant genes based on log2 fold change and adjusted p-value |
| `ma_plot.png`                | MA plot showing gene-wise fold changes vs mean expression |
| `heatmap_top_genes.png`      | Heatmap of the top 30 differentially expressed genes |
| `top_gene_count_plot.png`    | Boxplot showing expression counts for the most significant gene |
| `log2fc_density.png`         | Density plot showing distribution of log2 fold changes |
| `padj_histogram.png`         | Histogram of adjusted p-values from DESeq2 results |
| `top_genes_barplot.png`      | Barplot of top 10 upregulated and downregulated genes |

These figures are auto-generated in the analysis script and provide visual confirmation of trends, significant genes, and expression differences between sample groups.

##  How These Were Generated

All plots were created using R with the `ggplot2`, `EnhancedVolcano`, and `pheatmap` packages, based on DESeq2 results.
