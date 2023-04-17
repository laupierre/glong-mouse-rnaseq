## 230328_A00558_0209_BHN75WDSX5

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)


anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*TRM_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"

#write.xlsx (a, "star_gene_raw_counts.xlsx", rowNames=F)


annot <- a
annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

torm <- c("gene_name", "gene_type", "mgi_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
a <- a[ ,-1]

# remove sample GL9
a <- a[ ,grep ("GL9", colnames (a), invert=TRUE)]

sampleTable <- data.frame (matrix (nrow=dim (a)[2], ncol=2))
colnames (sampleTable) <- c("sample", "condition")
sampleTable$sample <- colnames (a)
row.names (sampleTable) <- colnames (a)
sampleTable$condition <- c("GL15", "GL15", "GL15", rep ("GL3", 6), rep ("GL15", 2), rep ("WT", 6))
sampleTable

idx <- match (sampleTable$sample, colnames (a))
sampleTable <- sampleTable[idx, ]
sampleTable$condition <- factor (sampleTable$condition)
stopifnot (sampleTable$sample == colnames (a))


dds <- DESeqDataSetFromMatrix(countData = round (a), colData = sampleTable, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds


## GL15 vs WT
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "GL15", "WT"))

res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Snca", ] 
write.xlsx (res, "GLONG_15monhs_vs_WT_2023.xlsx", rowNames=F)



## GL3 vs WT

res <- results(dds, contrast=c("condition", "GL3", "WT"))

res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Snca", ] 

write.xlsx (res, "GLONG_3months_vs_WT_2023.xlsx", rowNames=F)


## PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "sample"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, label=sample)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  		geom_text_repel()  + 
		coord_fixed () 

ggsave ("PCA plot.pdf")



## Heatmap plot

df <- as.data.frame(colData(dds)[,c("condition","sample")])
df$sample <- gsub ("_SND", "", df$sample)

# Th, DDC, NR4A2, SLC6A3
select <- c("ENSMUSG00000000214.12", "ENSMUSG00000020182.17", "ENSMUSG00000026826.14", "ENSMUSG00000021609.7", "ENSMUSG00000025889.14")

pdf ("Heatmap plot.pdf")
pheatmap(log2 (counts(dds,normalized=TRUE)+1) [row.names (counts(dds)) %in% select,], cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, annotation_col=df)
dev.off ()


## GL15 vs GL3

res <- results(dds, contrast=c("condition", "GL15", "GL3"))

res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

# Sanity check
res[res$gene_name == "Snca", ] 

write.xlsx (res, "GLONG_15months_vs_GLONG_3months_2023.xlsx", rowNames=F)












