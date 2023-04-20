library (openxlsx)
library (gprofiler2)

## The enrichment analysis takes only significantly differentially expressed genes

# 3 months GLONG mouse
res <- read.xlsx ("GLONG_3months_vs_WT_2023.xlsx")

res <- res[res$padj <= 0.05, ]
genes.up <- gsub ("\\..*", "", res$Geneid[res$log2FoldChange >0]
genes.down <- res$Geneid[res$log2FoldChange <0]

gp_up <- gost(query = genes.up, organism = "mmusculus", significant= FALSE, numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE)
# remove the evcodes and the intersections
gp_up <- gp_up$result[ ,-15]
gp_up <- gp_up[ ,-15]
head(gp_up)

