BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("RColorBrewer")
getwd()
setwd("C:/Users/msc2/Desktop/Assignments/TP/Prac_5_Deseq2")
count <- as.matrix(read.csv("count_data.csv",sep=",",row.names="gene_id"))
head(count)
colnames(count)
meta <- read.csv("metadata.csv", row.names=1)
head(meta)
rownames(meta)
colnames(meta)
all(rownames(meta) == colnames(count))
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = count, colData =meta, design = ~condition)
dim(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dim(dds)
dds <- DESeq(dds)
results <- results(dds)
results <- results[order(results$padj),]
head(results)
write.csv(res,"Treated_vs_control_DEG.csv")
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix)
colnames(sampleDistMatrix)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
plotPCA(vsd, intgroup="condition")
res <- results(dds, contrast=c("condition","treated","untreated"), alpha=0.05) 
summary(res)
pdf("MA_plot.pdf", width=6, height=6)
plotMA(res, ylim=c(-2,2))
dev.off()
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
results <- as.data.frame(res)
results <- na.omit(results)
results$significant <- ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 1,
                              "Significant", "Not Significant")

ggplot(results, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal()
