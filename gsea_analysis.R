library(msigdbr)
library (openxlsx)
library (fgsea)

# BP from GO
pathways <- msigdbr("mouse", category="C5", subcategory = "GO:BP")
pathways <- split (as.character (pathways$ensembl_gene), pathways$gs_name)

# 3 months
res <- read.xlsx ("GLONG_3months_vs_WT_2023.xlsx")

res <- res[order (res$log2FoldChange), ]
ranks <- res$log2FoldChange
names (ranks) <- gsub ("\\..*", "", res$Geneid)
  
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = ranks,
                  minSize  = 10,
                  maxSize  = 500)
  
fgseaRes <- fgseaRes[ ,c(1:6)]
fgseaRes <- data.frame (fgseaRes[order (fgseaRes$padj), ])
fgseaRes <- fgseaRes[order (fgseaRes$padj), ]
fgseaRes <- fgseaRes[ ,-4]
fgseaRes <- fgseaRes[ ,-4]
head (fgseaRes)

# In 8, we define the leading-edge subset to be those genes in the gene set S  
# that appear in the ranked list L at, or before, the point where the running sum reaches its maximum deviation from zero.

write.xlsx (fgseaRes, "gsea_pathways_3months_glong_vs_WT.xlsx", sep="") , rowNames=F)

