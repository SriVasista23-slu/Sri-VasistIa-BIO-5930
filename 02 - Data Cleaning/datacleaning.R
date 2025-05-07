# Load gene expression count matrix
counts <- read.table("GSE81089_readcounts_featurecounts1.tsv", 
                     header = TRUE, sep = "\t", row.names = 1)

# Load metadata and parse condition info
meta_lines <- readLines("GSE81089_series_matrix (2).txt")
char_line <- grep("!Sample_characteristics_ch1", meta_lines, value = TRUE)[1]
raw_conditions <- unlist(strsplit(char_line, "\t"))[-1]
group <- ifelse(grepl("T\\\"?$", raw_conditions), "Tumor", "Normal")

# Create metadata DataFrame
sample_info <- data.frame(
  sample = colnames(counts),
  condition = factor(group)
)

# Match metadata sample IDs with count matrix
stopifnot(all(colnames(counts) == sample_info$sample))

# Remove genes with missing values
counts_clean <- counts[complete.cases(counts), ]

# Filter genes with total counts > 10
counts_clean <- counts_clean[rowSums(counts_clean) > 10, ]

# Save cleaned count matrix
write.csv(counts_clean, "data/cleaned_counts.csv", row.names = TRUE)
write.csv(sample_info, "data/cleaned_metadata.csv", row.names = FALSE)
