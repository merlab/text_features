# library(writexl)
library(caret)
library(openxlsx)
set.seed(123)
cells <- readRDS("./cells.rds")
inDir <- "./result/ml_model_list_pathways/"
files <- list.files(inDir, pattern = ".*-.*.rds")
files <- gsub("\\.rds", "", files)
table_rf_test <- data.frame()
table_rf_train <- data.frame()
table_glm_test <- data.frame()
table_glm_train <- data.frame()
wb <- createWorkbook()
r <- c(good = 0, bad = 0)
g <- c(good = 0, bad = 0)


calculateMetrics <- function(pred, obs) {
  df <- data.frame(pred = pred, obs = obs)
  df <- df[!is.na(df$pred), ]
  df <- df[!is.na(df$obs), ]
  pred <- df$pred
  obs <- df$obs
  out <- c()
  out$pearsonCor <- cor(pred, obs, method = "pearson")
  out$spearmanCor <- cor(pred, obs, method = "spearman")
  out$kendallCor <- cor(pred, obs, method = "kendall")
  predNorm <- (pred - mean(pred)) / sd(pred)
  obsNorm <- (obs - mean(obs)) / sd(obs)
  out$RMSE <- RMSE(pred = predNorm, obs = obsNorm)
  out$MSE <- mean((obsNorm - predNorm)^2)
  out$MAE <- MAE(pred = predNorm, obs = obsNorm)
  return(out)
}

trainIn <- readRDS("./result/ml_model_list/ccle_ctrpv2_rnaseq.rds")
trainOut <- readRDS("./result/ml_model_list/ccle_ctrpv2_aac.rds")
testIn <- readRDS("./result/ml_model_list/gdse_rnaseq.rds")
testOut <- readRDS("./result/ml_model_list/gdse_aac.rds")

for (i in files) {
  mlModel <- readRDS(sprintf("%s/%s.rds", inDir, i))
  name <- gsub("-rf|-glmnet", "", i)
  drug <- strsplit(name, "-")[[1]][1]
  genes <- mlModel$genes

  trainPred <- predict(mlModel, na.omit(trainIn[, genes]))
  metrics <- calculateMetrics(pred = trainPred, obs = trainOut[names(trainPred), drug])
  vtrain <- unlist(c(drug, "pathway", metrics))

  testPred <- predict(mlModel, na.omit(testIn[, genes]))
  metrics <- calculateMetrics(pred = testPred, obs = testOut[names(testPred), drug])
  vtest <- unlist(c(drug, "pathway", metrics))

  if (grepl("-rf", i)) {
    table_rf_test <- rbind(vtest, table_rf_test)
    table_rf_train <- rbind(vtrain, table_rf_train)
  } else {
    table_glm_test <- rbind(vtest, table_glm_test)
    table_glm_train <- rbind(vtrain, table_glm_train)
  }

  v <- c()
}
coln <- c(
  "drug", "feature selection method",
  "pearsonCor",
  "spearmanCor",
  "kendallCor",
  "RMSE",
  "MSE",
  "MAE"
)
l <- list(
  "ElasticNet perf metrics-train" = table_glm_train,
  "RandomForest perf metrics-train" = table_rf_train,
  "ElasticNet perf metrics-test" = table_glm_test,
  "RandomForest perf metrics-test" = table_rf_test
)
l <- lapply(l, function(x) {
  colnames(x) <- coln
  for (i in 3:ncol(x)) {
    x[, i] <- as.numeric(x[, i])
  }
  # x <- x[order(abs(x$distanceToTarget), decreasing = FALSE), ]
  return(x)
})



warm1Style <- createStyle(fontColour = "#FFFFFF", bgFill = "#FF0000")
wb <- createWorkbook()
for (i in names(l)) {
  x <- l[[i]]
  addWorksheet(wb, i, gridLines = TRUE)
  writeData(wb, i, x)
}
saveWorkbook(wb, "./mlModelMetrics_pathway.xlsx", overwrite = TRUE)
print("done")
