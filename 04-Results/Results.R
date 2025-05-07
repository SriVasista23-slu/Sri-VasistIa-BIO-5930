# Load required packages
library(DESeq2)

# Assume dds and res already exist from the analysis
# dds: DESeqDataSet object
# res: DESeq2 results object

# Clean and sort results
res_clean <- na.omit(res)
res_sorted <- res_clean[order(res_clean$log2FoldChange), ]

# Save all significant DEGs to CSV
sig_genes <- res_clean[res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1, ]
write.csv(sig_genes, file = "results/significant_DEGs_GSE81089.csv", row.names = TRUE)

# Select top 10 upregulated and top 10 downregulated genes
top_down <- head(res_sorted, 10)
top_up <- tail(res_sorted, 10)
top_genes <- rbind(top_down, top_up)
write.csv(top_genes, file = "results/top_genes_log2FC.csv", row.names = TRUE)

# Optional: Save a summary of DEG counts
up_count <- nrow(sig_genes[sig_genes$log2FoldChange > 1, ])
down_count <- nrow(sig_genes[sig_genes$log2FoldChange < -1, ])
summary_text <- paste("Upregulated genes:", up_count,
                      "\nDownregulated genes:", down_count)
write(summary_text, file = "results/summary_DEG_counts.txt")
