library (openxlsx)

wt <- read.xlsx ("/Volumes/king/glong2023/2023/230328_A00558_0209_BHN75WDSX5_glong_wt/Comparison of GLONG vs WT at 3 and 15 months.xlsx")
ola <- read.xlsx ("/Volumes/king/glong2023/2022/220120_A00558_0160_AH5V77DSX3_glong_ola/Comparison of GLONG vs OLA at 3, 6 and 15 months.xlsx")

mco <- merge (ola, wt, by="Geneid", all.x=TRUE, all.y=TRUE)
mco <- mco [order (mco$`padj-3months-2022`), ]
colnames (mco)[colnames (mco) == "stats"] <- "stats-GLONGvsOLA"
colnames (mco)[colnames (mco) == "trend"] <- "trend-GLONGvsOLA"

write.xlsx (mco, "Comparison of GLONG vs OLA vs WT at each time point.xlsx", rowNames=F)