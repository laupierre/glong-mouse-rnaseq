library(msigdbr)
library (openxlsx)
library (fgsea)

# The Gene Set Enrichment Analysis (GSEA) differs from Gene Ontology enrichment analysis in that it considers all genes  
# in contrast to taking only significantly differentially expressed genes. 



mygsea <- function (res, outfile) {

## BP from GO
#pathways <- msigdbr("mouse", category="C5", subcategory = "GO:BP")
#pathways <- split (as.character (pathways$ensembl_gene), pathways$gs_name)

# Reactome
pathways <- msigdbr("mouse", category="C2", subcategory = "CP:REACTOME")
pathways <- split (as.character (pathways$ensembl_gene), pathways$gs_name)

# ranks are from lowest to highest
res <- res[order (res$log2FoldChange), ]
ranks <- res$log2FoldChange
names (ranks) <- gsub ("\\..*", "", res$Geneid)
  
fgseaRes <- fgsea(pathways = pathways, 
                  stats    = ranks,
                  minSize  = 10,
                  maxSize  = 500)
  
# see 20 pathways ordered by pval  
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam=0.5)


# Enrichment plot of a specific top pathway (up and down)
library (ggplot2)
plotEnrichment(pathways[[topPathwaysUp[[1]] ]], ranks) + ggtitle( topPathwaysUp[[1]] )

plotEnrichment(pathways[[topPathwaysDown[[1]] ]], ranks) + ggtitle( topPathwaysDown[[1]] )


# make leading edge more human-readable
library (org.Mm.eg.db)

fgseaRes[ , leadingEdge := mapIdsList(x=org.Mm.eg.db, 
                                      keys=leadingEdge,
                                      keytype="ENSEMBL", 
                                      column="SYMBOL")]


# save fgsea results in a file
fgseaRes <- fgseaRes[ ,c(1:6,8)]
fgseaRes <- data.frame (fgseaRes[order (fgseaRes$padj), ])
fgseaRes <- fgseaRes[order (fgseaRes$padj), ]
fgseaRes <- fgseaRes[ ,-4]
fgseaRes <- fgseaRes[ ,-4]
head (fgseaRes)

# In 8, we define the leading-edge subset to be those genes in the gene set S  
# that appear in the ranked list L at, or before, the point where the running sum reaches its maximum deviation from zero.

# Positive NES usually correspond to up-regulated gene sets and negative NES to down-regulated ones

write.xlsx (fgseaRes, outfile, rowNames=F)
}


# 3 months GLONG mouse
res <- read.xlsx ("GLONG_3months_vs_WT_2023.xlsx")
mygsea (res, "gsea_reactome_3months_glong_vs_WT.xlsx")







