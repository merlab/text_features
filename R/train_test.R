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

trainIndex <- createDataPartition(y, p = .8, 
                                  list = TRUE, 
                                  times = 3)

tgrid <- expand.grid(alpha=seq(0, 1, 0.2),
                     lambda=seq(0, 10, 1))

tcontrol <- trainControl(method="repeatedcv",
                        number=10,
                        repeats = 15,
                        search="grid",
                        savePredictions ="all",
                        allowParallel = TRUE,
                        verboseIter=TRUE)

pred_sample <- data.frame()

for(res in names(trainIndex))
{
  trIndx <- trainIndex[[res]]
  tsIndx <- setdiff(1:nrow(x), trIndx)
  
  train_result_sample <- train(x=x[trIndx, ], y=y[trIndx],
                               method="glmnet",
                               #preProcess=c("center", "scale"),
                               maximize = TRUE,
                               tuneGrid=tgrid,
                               trControl=tcontrol)
  train_result_sample$trainingData <- NULL
  pred_sample_n <- data.frame(index = tsIndx,
                              predict=predict(train_result_sample, x[tsIndx, ]),
                              original=y[tsIndx], 
                              resample=res)
  pred_sample <- rbind(pred_sample, pred_sample_n)
}

saveRDS(train_result_sample, "model.rds")

bwplot(train_result_sample$resample$RMSE,  xlab="RMSE", ylab = "Erlotinib", main = "Model Accuracy by Drug")

