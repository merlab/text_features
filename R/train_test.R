library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
source("./summarizeData.R")
GDSC2 <- readRDS("../data/GDSC2.rds")


df = summarizeData(pSet=GDSC2, mDataType = "Kallisto_0.46.1.rnaseq", drug = "Erlotinib", sensitivity.measure="aac_recomputed")


df = df[1:200]
df = df[, !is.na(colData(df)$aac_recomputed)]

NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
df <- df[, !NAsamples]

gene_vars <- apply(assay(df), 1, var)
df <- df[gene_vars>0, ]

x <- t(assay(df))
y <- colData(df)$aac_recomputed

train_result_sample <- train(x=x,
                      y=y,
                      method="glmnet",
                      preProcess=c("center", "scale"),
                      maximize = TRUE,
                      tuneGrid=expand.grid(alpha=seq(0, 1, 0.2),
                                           lambda=seq(0, 10, 1)),
                      trControl=trainControl(method="repeatedcv",
                                             number=10,
                                             repeats = 15,
                                             search="grid",
                                             savePredictions ="all",
                                             allowParallel = TRUE,
                                             verboseIter=TRUE))
train_result_sample$trainingData <- NULL
saveRDS(train_result_sample, "model.rds")

pred_sample <- data.frame(predict=predict(train_result_sample, x), original=y)
bwplot(train_result_sample$resample$RMSE,  xlab="RMSE", ylab = "Erlotinib", main = "Model Accuracy by Drug")

