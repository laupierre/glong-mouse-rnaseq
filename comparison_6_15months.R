library (openxlsx)

m3 <- read.xlsx ("GLONG_3months_vs_WT_2023.xlsx")
colnames (m3) <- paste (colnames (m3), "-3months-2023", sep="")

m15 <- read.xlsx ("GLONG_15months_vs_WT_2023.xlsx")
colnames (m15) <- paste (colnames (m15), "-15months-2023", sep="")

mco <- merge (m3, m15, by.x="Geneid-3months-2023", by.y="Geneid-15months-2023")

mco <- mco[ ,grep ("type|external|description|lfcSE|stat|pvalue|mgi|SND" , colnames (mco), invert=TRUE)]
mco <- mco [order (mco$`padj-3months-2023`), ]
colnames (mco)[colnames (mco) == "Geneid-3months-2023"] <- "Geneid"

mco$stats <- ifelse (mco$`padj-3months-2023` < 0.05 & mco$`padj-15months-2023` < 0.05, "Yes","No")
mco$trend <- ifelse (mco$`log2FoldChange-3months-2023` < 0 & mco$`log2FoldChange-15months-2023` < 0, "Down","No")
mco$trend [mco$`log2FoldChange-3months-2023` > 0 & mco$`log2FoldChange-15months-2023` > 0] <- "Up"
write.xlsx (mco, "Comparison of GLONG vs WT at 3 and 15 months.xlsx", rowNames=F)
