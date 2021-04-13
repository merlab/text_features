library(caret)
d5 <- readRDS("../data/Bortezomib_500.rds")
#d1 <- readRDS("../data/Bortezomib_100.rds")

x = d5$X; y=d5$Y$aac

method = "glmnet"
tgrid <- expand.grid(alpha=seq(0, 1, 0.2), lambda=seq(0, 10, 1))
tcontrol <- trainControl(method="repeatedcv",
                         number= 4,
                         repeats = 4,
                         search="grid",
                         savePredictions ="all",
                         allowParallel = TRUE,
                         verboseIter=FALSE)

trainIndex <- createDataPartition(y, p = .8, list = TRUE, times = 5)
mod100List <- list(); mod500List <- list();

for(res in names(trainIndex))
{
  trIndx <- trainIndex[[res]]
  tsIndx <- setdiff(1:nrow(x), trIndx)

  preProcValues <- preProcess(x[trIndx, ], method = c("center", "scale"))
  trainTransformed <- predict(preProcValues, x[trIndx, ])
  testTransformed <- predict(preProcValues, x[tsIndx, ])

  gene_vars <- abs(apply(x[trIndx,], 2, var))
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  var100 <- names(gene_vars)[1:100]
  var500 <- names(gene_vars)[1:500]

  tr100 <- trainTransformed[, var100]
  ts100 <- testTransformed[ , var100]

  tr500 <- trainTransformed[, var500]
  ts500 <- testTransformed[ , var500]

  ##-----------
  getTrainMod <- function(trData, y, trIndx, tgrid, tcontrol, tsIndx, testData,res)
  {
    model <- train(x=trData, y=y[trIndx],method="glmnet", tuneGrid=tgrid,
                      trControl=tcontrol)
    model$trainingData <- NULL

    predx <- data.frame(index = tsIndx, pred=predict(model, testData),
                        obs=y[tsIndx], resample=res)
    performance <- postResample(pred = predx$pred, obs = predx$obs)
    performance['cor'] <- cor(predx$pred, predx$obs)

    return(list(model=model, pred=predx, performance=performance))
  }

  model100 <- getTrainMod(tr100, y, trIndx,tgrid, tcontrol, tsIndx, ts100, res)
  model500 <- getTrainMod(tr500, y, trIndx,tgrid, tcontrol, tsIndx, ts500, res)

  mod100List[[res]] <- model100
  mod500List[[res]] <- model500
}

mean(sapply(mod100List, function(i)i$performance['cor']))
## [1] 0.3469093

mean(sapply(mod500List, function(i)i$performance['cor']))
## [1] 0.3530966


