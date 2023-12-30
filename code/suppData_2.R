library(writexl)
# loading the text mining data
dir <- "./data/drug_text-features"
l <- list()
for (i in list.files(dir, patter = "*.rds")) {
    drug <- gsub("\\.rds", "", i)
    x <- readRDS(sprintf("%s/%s", dir, i))
    x <- x[order(x$FDR, decreasing = FALSE), ]
    l[[drug]] <- x
}
write_xlsx(l, "./result/Supplementary_Data_2_txtMiningGeneList.xlsx")
