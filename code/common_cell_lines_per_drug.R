source("./code/helper/summarizeData.R")
source("./code/helper/computeInteractionMatrix.R")
source("./code/helper/mRMR.R")
source("./code/helper/summarizeData_pharmacoGx.R")
library(PharmacoGx)
library(Xeva)
library(writexl)
removeCombTherapies <- function(v) {
  v <- v[!grepl(":", v)]
  v <- v[!grepl("^ML", v)]
  v <- v[!grepl("^MK", v)]
  v <- v[!grepl("^VER", v)]
  v <- v[!grepl("^VU", v)]
  v <- v[!grepl("^SKI", v)]
  v <- v[!grepl("^SJ", v)]
  v <- v[!grepl("[0-9]", v)]
  v <- v[!grepl("^N-", v)]
  v <- v[!grepl("^\\(-\\)-", v)]
  v <- v[!grepl("-leu-", v)]
  return(v)
}

getCellLines <- function(mat, drug) {
  eset <- summarizeDataPGx(
    mat,
    mDataType = "Kallisto_0.46.1.rnaseq",
    drug = drug
  )
  cData <- colData(eset)
  cells <- rownames(cData[!is.na(cData$aac_recomputed), ])
  return(cells)
}
ccle <- readRDS("./data/ccle_ctrpv2_aac.rds")
gdse <- readRDS("./data/gdse_aac.rds")

gdse <- readRDS("./data/PSet_GDSC2020.rds")
gdsedrugs <- rownames(drugInfo(gdse))
clinfo <- cellInfo(gdse)

ccle <- readRDS("./data/CCLE-CTRPv2_Kallisto_0.46.1.rnaseq.rds")
ccledrugs <- rownames(drugInfo(ccle))
clinfo <- cellInfo(ccle)
alldrugs <- unique(c(ccledrugs, gdsedrugs))
out <- data.frame()
f <- "./tmp.rds"
if (file.exists(f)) {
  out <- readRDS(f)
  alldrugs <- alldrugs[!alldrugs %in% out$drugName]
}
for (drug in alldrugs) {
  print(drug)
  if (drug %in% ccledrugs) {
    ccleCells <- getCellLines(ccle, drug)
  } else {
    ccleCells <- c()
  }
  if (drug %in% gdsedrugs) {
    gdseCells <- getCellLines(gdse, drug)
  } else {
    gdseCells <- c()
  }
  res <- c(
    drug, length(ccleCells),
    length(gdseCells),
    length(intersect(ccleCells, gdseCells))
  )
  out <- rbind(res, out)
  colnames(out) <- c("drugName", "n_cells_CCLE", "n_cells_GDSE", "n_common_cells")
  saveRDS(out, f)
}
colnames(out) <- c("drugName", "n_cells_CCLE", "n_cells_GDSE", "n_common_cells")
out$"n_cells_CCLE" <- as.numeric(out$"n_cells_CCLE")
out$"n_cells_GDSE" <- as.numeric(out$"n_cells_GDSE")
out$"n_common_cells" <- as.numeric(out$"n_common_cells")
l <- list("allDrugs" = out, "monotherapies" = out[out$drugName %in% removeCombTherapies(out$drugName), ])
# alldrugs <-
write_xlsx(l, "./result/common_cell_lines_per_drug.xlsx")
print("done")
print("done")
