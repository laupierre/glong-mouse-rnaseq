library (openxlsx)

m3 <- read.xlsx ("GLONG_3months_vs_WT_2023.xlsx")
colnames (m3) <- paste (colnames (m3), "-3months-WT-2023", sep="")
#reverse the logFC
m3$`log2FoldChange-3months-WT-2023` <- -1 * m3$`log2FoldChange-3months-WT-2023`

m15 <- read.xlsx ("GLONG_15months_vs_WT_2023.xlsx")
colnames (m15) <- paste (colnames (m15), "-15months-WT-2023", sep="")
#reverse the logFC
m15$`log2FoldChange-15months-WT-2023` <- -1 * m15$`log2FoldChange-15months-WT-2023`


mco <- merge (m3, m15, by.x="Geneid-3months-WT-2023", by.y="Geneid-15months-WT-2023")

mco <- mco[ ,grep ("type|external|description|lfcSE|stat|pvalue|mgi|SND" , colnames (mco), invert=TRUE)]
mco <- mco [order (mco$`padj-3months-WT-2023`), ]
colnames (mco)[colnames (mco) == "Geneid-3months-WT-2023"] <- "Geneid"

mco$stats <- ifelse (mco$`padj-3months-WT-2023` < 0.05 & mco$`padj-15months-WT-2023` < 0.05, "Yes","No")
mco$trend <- ifelse (mco$`log2FoldChange-3months-WT-2023` < 0 & mco$`log2FoldChange-15months-WT-2023` < 0, "Down","No")
mco$trend [mco$`log2FoldChange-3months-WT-2023` > 0 & mco$`log2FoldChange-15months-WT-2023` > 0] <- "Up"
colnames (mco)[colnames (mco) == "stats"] <- "stats-WTvsGLONG"
colnames (mco)[colnames (mco) == "trend"] <- "trend-WTvsGLONG"

write.xlsx (mco, "Comparison of GLONG vs WT at 3 and 15 months.xlsx", rowNames=F)
