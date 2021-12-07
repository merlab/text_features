library(caret)
#d5 <- readRDS("../data/Bortezomib_500.rds")
#d1 <- readRDS("../data/Bortezomib_100.rds")

d5 <- readRDS("~/Downloads/Cediranib_500.rds")
x = d5$X; y=d5$Y$aac

method = "glmnet"
tgrid <- expand.grid(alpha=seq(0.0, 1, 0.1),
                     lambda = 10^seq(-3, 3, length = 10))

tcontrol <- trainControl(method = "cv", #"repeatedcv",
                         number= 4, #repeats = 4,
                         search="grid", allowParallel = TRUE,
                         verboseIter=FALSE,
                         savePredictions=c("all", "final", "none")[3],
                         returnResamp   =c("all", "final", "none")[3],
                         returnData = FALSE)

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
  train_elastic_net_model <- function(trData, y, trIndx, tgrid, tcontrol, tsIndx, testData,res)
  {
    model <- train(x=trData, y=y[trIndx],method="glmnet", trControl=tcontrol,
                   tuneGrid=tgrid #tuneLength = 10 #
                   )
    #model$trainingData <- NULL ##no need as returnData=F
    model$finalModel$call <- NULL
    ##----feature names are in model$finalModel$xNames ------

    predx <- data.frame(index = tsIndx, pred=predict(model, testData),
                        obs=y[tsIndx], resample=res)
    performance <- postResample(pred = predx$pred, obs = predx$obs)
    performance['cor'] <- cor(predx$pred, predx$obs)

    return(list(model=model, pred=predx, performance=performance))
  }

  model100 <- train_elastic_net_model(tr100, y, trIndx,tgrid, tcontrol, tsIndx, ts100, res)
  model500 <- train_elastic_net_model(tr500, y, trIndx,tgrid, tcontrol, tsIndx, ts500, res)

  mod100List[[res]] <- model100
  mod500List[[res]] <- model500
}

mean(sapply(mod100List, function(i)i$performance['cor']))
## [1] 0.4272153

mean(sapply(mod500List, function(i)i$performance['cor']))
## [1] 0.4641324


