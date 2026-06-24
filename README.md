# RNA-Seq Differential Gene Expression Analysis Pipeline

This repository contains an end-to-end workflow for analyzing RNA-seq raw count data using the **DESeq2** Bioconductor package. The dataset is derived from a study by Bellat et al. (2020) examining the transcriptional profile of MDA-MB-468 triple-negative breast cancer cells treated with **salinomycin** versus a PBS control.

---

## 📊 Dataset Overview


**Organism:** *Homo sapiens* (Human) 

**Platform:** Illumina HiSeq 4000 

**Reference Genome:** GRCh37 mapped via STAR 

**Experimental Design:** 3 Replicates of Untreated Controls vs. 3 Replicates of Salinomycin-Treated samples (24-hour time point).

**Data Accession:** Raw count data is archived under GEO accession [GSE135514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135514).



---

## 🛠️ Typical Analysis Workflow

```
   [ FASTQ ] ──> Quality Control (FastQC)
       │
       ▼
 [ Pre-processing ] ──> Adapter/Quality Trimming (Trimmomatic)
       │
       ▼
   [ Alignment ] ──> Mapping to Reference Genome (STAR/TopHat)
       │
       ▼
 [ Quantitation ] ──> Count Generation per Gene (HTSeq/Cufflinks)
       │
       ▼
 [ Downstream ] ──> Normalization & Differential Expression (DESeq2)

```

---

## 🚀 Getting Started

### Prerequisites & Installation

Ensure you have **R** installed. Run the following block to install the required Bioconductor and CRAN packages:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "RColorBrewer", "EnhancedVolcano"))
install.packages("ggplot2")

```

### Input Files Required

1. 
`count_data.csv`: Matrix featuring raw integers of reads mapped to individual gene Ensembl IDs.


2. 
`metadata.csv`: Annotation mapping sample titles to their experimental `condition` group.



---

## 💻 Core Implementation Script

```R
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("EnhancedVolcano")

# 1. Load Data
[cite_start]count <- as.matrix(read.csv("count_data.csv", sep=",", row.names="gene_id")) [cite: 163]
[cite_start]meta <- read.csv("metadata.csv", row.names=1) [cite: 168]

# Verify alignment of row names and column names
[cite_start]if(!all(rownames(meta) == colnames(count))) stop("Metadata and Counts columns do not match!") [cite: 172, 173]

# 2. Setup DESeqDataSet object
[cite_start]dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~condition) [cite: 176]

# 3. Pre-filtering low count genes
[cite_start]keep <- rowSums(counts(dds)) >= 10 [cite: 179]
[cite_start]dds <- dds[keep,] [cite: 180]

# 4. Differential Expression Analysis
[cite_start]dds <- DESeq(dds) [cite: 182]
[cite_start]res <- results(dds, contrast=c("condition", "treated", "untreated"), alpha=0.05) [cite: 207]
[cite_start]res <- res[order(res$padj),] [cite: 184]

# Save outputs
[cite_start]write.csv(as.data.frame(res), "Treated_vs_control_DEG.csv") [cite: 186]

# 5. Data Transformation for Visualizations
[cite_start]vsd <- vst(dds, blind=FALSE) [cite: 188]

```

---

## 📈 Quality Control & Exploratory Data Analysis

### Sample Distance Clustering (Heatmap)

```R
[cite_start]sampleDists <- dist(t(assay(vsd))) [cite: 195]
[cite_start]sampleDistMatrix <- as.matrix(sampleDists) [cite: 196]
[cite_start]colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) [cite: 199]
[cite_start]pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors) [cite: 200, 201, 202, 203]

```

* 
**Result Interpretation:** The three control samples cluster tightly together, separate from the salinomycin-treated group. This establishes solid experimental reproducibility and proves the drug response causes robust transcriptional remodeling rather than technical noise.



### Principal Component Analysis (PCA)

```R
[cite_start]plotPCA(vsd, intgroup="condition") [cite: 204, 205]

```

* 
**Result Interpretation:** **PC1 accounts for 96% of total variance**. The clear separation along the primary component confirms that the drug treatment condition is the dominant source of variation across the sample set.



---

## 📋 Statistical Results Summary

Out of **19,597 genes** evaluated with non-zero expression values:

| Category | Count | Percentage | Details |
| --- | --- | --- | --- |
| **Upregulated ($LFC > 0$)** | 2,369 

 | ~12% 

 | Significantly higher in Salinomycin group ($p_{adj} < 0.05$) 

 |
| **Downregulated ($LFC < 0$)** | 2,222 

 | ~11% 

 | Significantly lower in Salinomycin group ($p_{adj} < 0.05$) 

 |
| **Low Counts Filtered** | 4,939 

 | ~25% 

 | Extracted due to mean normalized count $< 15$ 

 |
| **Outliers** | 0 

 | 0% 

 | None flagged by Cook's distance filter 

 |

The balanced distribution (~4,600 overall differentially expressed genes) reveals massive, systematic bidirectional genetic shifts triggered by the treatment.

---

## 🎨 Downstream Visualizations

### MA Plot

```R
[cite_start]plotMA(res, ylim=c(-2,2)) [cite: 209, 211]

```

Plots the relationship between mean normalized expression level ($X$-axis) against $\log_2$ fold changes ($Y$-axis). Significantly shifting elements are highlighted in blue, exhibiting higher dispersion at low expression thresholds.

### Volcano Plot

```R
[cite_start]results_df <- na.omit(as.data.frame(res)) [cite: 215, 216]
[cite_start]results_df$significant <- ifelse(results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1, "Significant", "Not Significant") [cite: 217, 218]

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  [cite_start]theme_minimal() [cite: 219, 220, 221, 222, 223]

```

Displays statistical significance ($-\log_{10} p_{adj}$) against magnitude of change ($\log_2\text{FoldChange}$) , immediately isolating candidate markers that meet both critical statistical and biological thresholds.
