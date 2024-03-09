# library(writexl)
library(caret)
library(openxlsx)
set.seed(123)
cells <- readRDS("./cells.rds")
commonCellLines <- intersect(cells[[1]], cells[[2]])
inDir <- "./result/ml_model_list/"
files <- list.files(inDir, pattern = ".*-.*.rds")
files <- gsub("\\.rds", "", files)
table_rf <- data.frame()
table_glm <- data.frame()
wb <- createWorkbook()


mlMethods <- c("var-100", "var-500", "cor-500", "text-mining", "L1000", "L1000-tm")

for (i in files) {
  mlModel <- readRDS(sprintf("%s/%s.rds", inDir, i))
  name <- gsub("-rf|-glmnet", "", i)
  drug <- strsplit(name, "-")[[1]][1]
  genes <- mlModel$genes
  geneFilterMethod <- paste0(unlist(strsplit(name, "-")[[1]][-1]), collapse = "-")
  tune <- mlModel$bestTune
  seed <- mlModel$seed
  type <- mlModel$modelType
  genes <- paste0(genes, collapse = ",")


  if (grepl("-rf-", i)) {
    idx <- nrow(table_rf) + 1
    table_rf[idx, "drug"] <- drug
    table_rf[idx, "feature selection method"] <- geneFilterMethod
    table_rf[idx, names(tune)] <- tune
    table_rf[idx, "seed"] <- seed
    table_rf[idx, "type"] <- type
    table_rf[idx, "selected genes"] <- genes
  } else {
    idx <- nrow(table_glm) + 1
    table_glm[idx, "drug"] <- drug
    table_glm[idx, "feature selection method"] <- geneFilterMethod
    table_glm[idx, names(tune)] <- tune
    table_glm[idx, "maximize"] <- mlModel$maximize
    table_glm[idx, "seed"] <- seed
    table_glm[idx, "type"] <- type
    table_glm[idx, "selected genes"] <- genes
    names(mlModel)
    str(mlModel)
    stop()
  }

  v <- c()
}

warm1Style <- createStyle(fontColour = "#FFFFFF", bgFill = "#FF0000")
wb <- createWorkbook()
addWorksheet(wb, "ElasticNet parameters", gridLines = TRUE)
writeData(wb, "ElasticNet parameters", table_glm)
addWorksheet(wb, "RandomForest parameters", gridLines = TRUE)
writeData(wb, "RandomForest parameters", table_rf)

saveWorkbook(wb, "./parameters.xlsx", overwrite = TRUE)
print("done")
