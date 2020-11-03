library(caret)
set.seed(3456)
x=iris[,1:4]; y=iris$Species
trainIndex <- createDataPartition(y, p = .8, 
                                  list = TRUE, 
                                  times = 3)
#head(trainIndex)
pred_sample <- data.frame()
for(res in names(trainIndex))
{
  trIndx <- trainIndex[[res]]
  tsIndx <- setdiff(1:nrow(x), trIndx)
  
  train_result_sample <- train(x=x[trIndx, ], y=y[trIndx],
                               method="glmnet",
                               #preProcess=c("center", "scale"),
                               maximize = TRUE,
                               tuneGrid=expand.grid(alpha=seq(0, 1, 0.2),
                                                    lambda=seq(0, 10, 1)),
                               trControl=trainControl(method="repeatedcv",
                                                      number=4,
                                                      search="grid",
                                                      savePredictions ="all",
                                                      allowParallel = TRUE,
                                                      verboseIter=TRUE))
  train_result_sample$trainingData <- NULL
  pred_sample_n <- data.frame(index = tsIndx,
                              predict=predict(train_result_sample, x[tsIndx, ]),
                              original=y[tsIndx], 
                              resample=res)
  pred_sample <- rbind(pred_sample, pred_sample_n)
}

saveRDS(pred_sample, "pred_result.rds")