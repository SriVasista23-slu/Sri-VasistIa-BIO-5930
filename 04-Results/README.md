# Results Folder

This folder contains the key output files generated from the differential expression analysis of lung cancer samples using the GSE81089 dataset. Each file represents a specific outcome from different stages of the statistical analysis pipeline.

##  1. `significant_DEGs_GSE81089.csv`
- **Description:** A comprehensive list of genes identified as significantly differentially expressed between tumor and normal lung tissue samples.
- **Criteria:** Genes were selected based on adjusted p-value < 0.05 and an absolute log2 fold change greater than 1.
- **Contents:** Includes gene IDs, base mean expression, log2 fold changes, standard errors, raw p-values, adjusted p-values, and other DESeq2-generated metrics.

##  2. `top_genes_log2FC.csv`
- **Description:** This file contains the top 10 genes with the highest increase and the top 10 with the largest decrease in expression (based on log2 fold change).
- **How it was generated:** Genes were ranked by log2 fold change after removing entries with missing values.
- **Use case:** These genes are featured in key plots (e.g., bar plots) and are candidates for downstream biological interpretation or validation.

##  3. `summary_DEG_counts.txt` **
- **Description:** A short summary text that shows how many genes were significantly upregulated or downregulated in tumor samples compared to normal tissue.
- **Purpose:** Helps provide a quick snapshot of the overall results, useful for summaries and figure legends.

---

All output files were derived from a rigorous RNA-Seq analysis pipeline built using DESeq2 in R.
