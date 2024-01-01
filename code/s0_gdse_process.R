# purpose: process datasets for CCLE into CSV files for deep learning
source("./R/helper/summarizeData.R")
source("./R/helper/computeInteractionMatrix.R")
source("./R/helper/mRMR.R")
source("./R/summarizeData_pharmacoGx.R")
library(PharmacoGx)
library(Xeva)
gdse <- readRDS("./data/PSet_GDSC2020.rds")
# gdse <- updateObject(gdse)
drinfo <- drugInfo(gdse)
drugs <- rownames(drinfo)
txtmining_drugs <- list.files("./data/drug_text-features/", pattern = "*.rds")
txtmining_drugs <- gsub("\\.rds", "", txtmining_drugs)
txtmining_drug <- txtmining_drugs[txtmining_drugs %in% drugs]
# clinfo <- cellInfo(ccle)

allRNA <- data.frame()
allAAC <- data.frame()
for (i in txtmining_drugs) {
  print(i)
  eset <- summarizeDataPGx(
    gdse,
    mDataType = "Kallisto_0.46.1.rnaseq",
    drug = i
  )
  print("eset loaded")
  # extract rnaseq data
  rData <- rowData(eset)
  rData <- rData[which(rData$"gene_type" == "protein_coding"), ]
  rData <- rData[!duplicated(rData$gene_name), ]
  cData <- colData(eset)
  eset <- eset[rownames(rData), ]
  rownames(eset) <- rData$gene_name
  rna <- as.data.frame(assays(eset)[["exprs"]])
  rna$gene <- rownames(rna)
  print("now merging rnaseq")
  if (nrow(allRNA) > 0) {
    rna <- rna[, !colnames(rna) %in% colnames(allRNA)]
    rna$gene <- rownames(rna)
    if (ncol(rna) > 0) {
      allRNA <- merge(rna, allRNA, by = "gene", all = TRUE)
    }
  } else {
    allRNA <- rna
  }


  #
  resp <- data.frame(cell_line_name = rownames(cData))
  resp[, i] <- cData$aac_recomputed
  print("now merging aac")
  if (nrow(allAAC) > 0) {
    allAAC <- merge(resp, allAAC, by = "cell_line_name", all = TRUE)
  } else {
    allAAC <- resp
  }
  print(dim(allRNA))
  print(dim(allAAC))
}
allRNA <- as.data.frame(allRNA)
rownames(allRNA) <- allRNA$gene
allRNA$gene <- NULL
allRNA <- t(allRNA)
rownames(allAAC) <- allAAC$cell_line_name
allAAC$cell_line_name <- NULL
# allRNA$cell_line_name <- rownames(allRNA)
write.csv(allAAC, "./result/gdse_aac.csv", quote = FALSE)
write.csv(allRNA, "./result/gdse_rnaseq_tpm.csv", quote = FALSE)
saveRDS(allAAC, "./result/gdse_aac.rds")
saveRDS(allRNA, "./result/gdse_rnaseq.rds")
print("done")
