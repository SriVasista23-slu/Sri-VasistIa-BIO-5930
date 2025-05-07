# Data Cleaning

This folder contains the preprocessing steps performed before the differential expression analysis.

## Contents
- `data_cleaning.R`: R script used to filter and prepare the raw RNA-Seq count data.

## What Was Done
- Removed any genes with missing (NA) values.
- Filtered out low-expression genes by excluding those with total counts ≤ 10 across all samples.
- Verified that the sample names in the count matrix match the sample information metadata.
- Prepared a cleaned count matrix used as input for DESeq2.

These steps helped ensure the analysis focused only on informative and reliable genes.
