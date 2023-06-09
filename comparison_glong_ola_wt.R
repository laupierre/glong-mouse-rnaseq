library (openxlsx)

wt <- read.xlsx ("/Volumes/king/glong2023/2023/230328_A00558_0209_BHN75WDSX5_glong_wt/Comparison of GLONG vs WT at 3 and 15 months.xlsx")
ola <- read.xlsx ("/Volumes/king/glong2023/2022/220120_A00558_0160_AH5V77DSX3_glong_ola/Comparison of GLONG vs OLA at 3, 6 and 15 months.xlsx")

mco <- merge (ola, wt, by="Geneid", all.x=TRUE, all.y=TRUE)
colnames (mco)[colnames (mco) == "stats"] <- "stats-GLONGvsOLA"
colnames (mco)[colnames (mco) == "trend"] <- "trend-GLONGvsOLA"

mco$overall_stats <- ifelse (mco$`stats-GLONGvsOLA` == "Yes" & mco$`stats-WTvsGLONG` == "Yes", "Yes", "No")
mco$overall_trend <- ifelse (mco$`trend-GLONGvsOLA` == "Up" & mco$`trend-WTvsGLONG` == "Up", "Up", "No")
mco$overall_trend [mco$`trend-GLONGvsOLA` == "Down" & mco$`trend-WTvsGLONG` == "Down"] <- "Down"

## Annotation
anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

mco <- merge (mco, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE)
mco <- mco [order (mco$`padj-3months-2022`), ]


write.xlsx (mco, "Comparison of GLONG vs OLA vs WT at each time point.xlsx", rowNames=F)
