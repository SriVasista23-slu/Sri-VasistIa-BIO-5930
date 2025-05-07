# RNA-Seq Analysis of Lung Tumor vs. Normal Tissue (GSE81089)

**Course Project – BIOL-5930: Data Analysis**  
**Institution:** Saint Louis University  
**Students:** Ifthekar Hussain & Sri Vasista Talagampala  

---

##  Project Overview

This project explores transcriptomic differences between lung tumor tissues and matched normal tissues using RNA-Seq data from the public dataset **GSE81089**. The objective is to identify significantly altered genes in lung tumors and gain biological insights into cancer progression.

---

##  Dataset Summary

- **Source:** NCBI Gene Expression Omnibus (GEO)  
- **Accession Number:** GSE81089  
- **Samples:** 199 tumor, 19 normal  
- **Organism:** *Homo sapiens*  
- **Technology:** RNA-Seq (processed using featureCounts)

---

##  Methods & Tools

Analysis was performed in **R** using Bioconductor packages. Main steps included:

- **Preprocessing:** Sample alignment and removal of low-count genes (sum < 10)
- **Normalization:** Median-of-ratios method (via DESeq2)
- **Statistical Testing:** Wald test with Benjamini-Hochberg correction
- **DEG Criteria:** Adjusted p-value < 0.05 and |log2 Fold Change| > 1

###  Key R Packages:
- `DESeq2` – Differential gene expression analysis  
- `ggplot2`, `pheatmap`, `EnhancedVolcano` – Visualization  
- `biomaRt` – Gene annotation and ID conversion  

---

##  Key Findings

- **Significant DEGs:** ~6,500 genes identified  
- **Top Upregulated:** ENSG00000234854 (~800-fold increase)  
- **Top Downregulated:** OLFM4 (ENSG00000102837; ~32-fold decrease)  
- **PCA Plot:** Strong separation between tumor and normal tissues  
- **Heatmap:** Top 30 DEGs clustered clearly by sample condition  
- **MA & Volcano Plots:** Highlighted strong differential expression patterns  

---

##  Biological Insights

The findings reveal substantial gene expression changes in lung cancer.  
Notable patterns include:
- **Suppression of OLFM4**, a potential protective gene in healthy lung epithelium
- **Upregulation of uncharacterized transcripts**, possibly linked to tumor growth

These results point to potential biomarkers and targets for future studies.

---

##  Repository Structure

| File / Folder                  | Description                                              |
|-------------------------------|----------------------------------------------------------|
| `finalprojectdataanalysis.Rmd`| R Markdown notebook with the complete analysis pipeline |
| `README.md`                   | Main documentation (this file)                          |
| `figures/`                    | Visual outputs (PCA, MA plot, heatmap, volcano plot)     |
| `results/`                    | DEG result tables and top gene summaries                |
| `scripts/`                    | Modular R scripts for pipeline stages                   |

---

##  Academic Context

This project was completed as part of **BIOL-5930: Data Analysis** at Saint Louis University.  
All data used in this study are publicly available and this analysis was conducted solely for educational purposes.

---

##  Authors

**Ifthekar Hussain**  
M.S. Bioinformatics, Saint Louis University  

**Sri Vasista Talagampala**  
M.S. Bioinformatics, Saint Louis University  
