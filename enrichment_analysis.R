library (openxlsx)
library (gprofiler2)
library (org.Mm.eg.db)

## The enrichment analysis takes only the enrichment of significantly differentially expressed genes


my_gprof <- function (res, outfile_up, outfile_down) {

# Take all expressed genes for the background
all_genes <- gsub ("\\..*", "", res$Geneid)

# keep only the significant genes
results_sig <- subset(res, padj < 0.05)
# get the significant up-regulated genes
up <- subset(results_sig, log2FoldChange > 0)
# get the significant down-regulated genes
down <- subset(results_sig, log2FoldChange < 0)


## Accounting for the order of genes in the enrichment analysis
## See https://f1000research.com/articles/9-709

# order genes by log2FC (from the most up-regulated to the least up-regulated)
up_ordered <- up[order(up$log2FoldChange, decreasing = TRUE),]
genes.up <- gsub ("\\..*", "", up_ordered$Geneid)

# order genes by log2FC (from the most down-regulated to the least down-regulated)
down_ordered <- down[order(down$log2FoldChange, decreasing = FALSE),]
genes.down <- gsub ("\\..*", "", down_ordered$Geneid)


## Up analysis
gp_up <- gost(query = genes.up, organism = "mmusculus", ordered_query = TRUE, significant= TRUE, numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), 
              exclude_iea=TRUE, evcodes = TRUE, correction_method = "gSCS", custom_bg = all_genes)
              
# Make more human readable
gp_up$result$genes_intersection <- ""

for (i in (1:dim (gp_up$result)[1])) {

k <- gp_up$result$intersection [i]
k <- unlist (strsplit (k, split=","))
a <- paste (sort (unique (mapIds(org.Mm.eg.db, keys=k, column=c("SYMBOL"), keytype="ENSEMBL"))), collapse=",")
gp_up$result$genes_intersection[i] <- a
}

# remove the evcodes and the intersections
gp_up <- gp_up$result[ ,-15]
gp_up <- gp_up[ ,-15]
# collect the results
gp_up <- gp_up[ ,c("term_id", "term_name", "p_value", "significant", "intersection_size", "genes_intersection")]
head(gp_up)
write.xlsx (gp_up, outfile_up, rowNames=F)



## Down analysis
gp_down <- gost(query = genes.down, organism = "mmusculus", ordered_query = TRUE, significant= TRUE, numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), 
              exclude_iea=TRUE, evcodes = TRUE, correction_method = "gSCS", custom_bg = all_genes)

# Make more human readable
gp_down$result$genes_intersection <- ""

for (i in (1:dim (gp_down$result)[1])) {

k <- gp_down$result$intersection [i]
k <- unlist (strsplit (k, split=","))
a <- paste (sort (unique (mapIds(org.Mm.eg.db, keys=k, column=c("SYMBOL"), keytype="ENSEMBL"))), collapse=",")
gp_down$result$genes_intersection[i] <- a
}

# remove the evcodes and the intersections
gp_down <- gp_down$result[ ,-15]
gp_down <- gp_down[ ,-15]
# collect the results
gp_down <- gp_down[ ,c("term_id", "term_name", "p_value", "significant", "intersection_size", "genes_intersection")]
head(gp_down)
write.xlsx (gp_down, outfile_down, rowNames=F)
}


# 3 months GLONG mouse
res <- read.xlsx ("GLONG_3months_vs_WT_2023.xlsx")
my_gprof (res, "GLONG_3months_vs_WT_2023_goa_up_genes.xlsx", "GLONG_3months_vs_WT_2023_goa_down_genes.xlsx")


# 15 months GLONG mouse
res <- read.xlsx ("GLONG_15months_vs_WT_2023.xlsx")
my_gprof (res, "GLONG_15months_vs_WT_2023_goa_up_genes.xlsx", "GLONG_15months_vs_WT_2023_goa_down_genes.xlsx")

# 3 vs 15 months GLONG mouse
res <- read.xlsx ("GLONG_15months_vs_GLONG_3months_2023.xlsx")
my_gprof (res, "GLONG_15months_vs_GLONG_3months_2023_goa_up_genes.xlsx", "GLONG_15months_vs_GLONG_3months_2023_goa_down_genes.xlsx")









