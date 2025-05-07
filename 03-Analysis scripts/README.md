# Analysis Scripts

This folder contains the main script used to analyze the cleaned RNA-Seq data for differential expression between lung tumor and normal samples.

## Main File

- **analysis.R**: This script performs:
  - Import of cleaned count and metadata files
  - Creation of DESeq2 dataset
  - Differential gene expression analysis
  - Visualization of results including PCA, MA plot, volcano plot, barplot, and heatmap

## Figures

All output plots are stored in the `figures/` subfolder:
- PCA plot
- MA plot
- Volcano plot
- Barplot of top DEGs
- Heatmap of top 30 DEGs
- Count plot for top gene

These figures help interpret the gene expression differences between tumor and normal tissue.
