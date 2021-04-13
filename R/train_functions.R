source("./summarizeData.R")
source("./computeInteractionMatrix.R")
source("./mRMR.R")

generate_df <- function(pSet, mDataType, drugname){
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drugname, sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]

  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]

  #remove samples with no variance
  gene_vars <- apply(assay(df), 1, var)
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  gene_vars <- gene_vars[gene_vars>0]
  df <- df[names(gene_vars), ]

  ##------code to remove features that have low variance ------
  ##atleast 15% should be non zero
  cutoff <- floor(ncol(df)*0.15)
  unique_count <- apply(assay(df), 1, function(v) length(unique(round(v, 4))))
  gn2take <- unique_count > cutoff
  df <- df[gn2take, ]
  ##------------------------------------------------------
  return(df)
}

subset_by_feat <- function(df, drug, textmining = NULL , subset_size = 0, cutoff_method = "waterfall", L1000 = FALSE, var_count = 0){
  minedgenes <- readRDS(sprintf("./genes/%s.rds",drug))
  #select mined genes
  if (L1000 == TRUE){
    dfout = df[rowData(df)$gene_name %in% L1000_gene_list, ]
  } else if (textmining == TRUE){
    dfout = df[rowData(df)$gene_name %in% minedgenes$Symbol, ]
  } else if (textmining == FALSE) {
    dfout= df[!(rowData(df)$gene_name %in% minedgenes$Symbol), ]
  }
  
  if (subset_size != 0){
    dfout = dfout[1:subset_size]
  }

  #produce x and y, turn y to discrete
  x <- t(assay(dfout))
  
  if (var_count > 0){
    gene_vars <- abs(apply(x, 2, var))
    gene_vars <- sort(gene_vars, decreasing = TRUE)
    gene_vars <- gene_vars[1:var_count]
    x <- x[, names(gene_vars)]
  }
  
  y <- colData(dfout)$aac_recomputed
  if (cutoff_method == "waterfall"){
    cutoff <- callingWaterfall(y, "AUC")
  } else {
    cutoff <- 0.2
  }
  y2 <- ifelse(y >= cutoff,"sensitive","resistant")

  #print class counts for y
  counts <- table(y2)
  print(sprintf("Resistant: %d", counts[1]))
  print(sprintf("Sensitive: %d", counts[2]))

  output <- list("X"=x, "Y" = list("class"=y2, "aac" = y))
  return(output)
}

#second function trains model based on x and y
trainmodel <- function(x,y,name,type, method, ft = -100, var_count = -100){
  set.seed(1)
  if(method == "glmnet"){
    tgrid <- expand.grid(alpha=seq(0, 1, 0.2),
                         lambda=seq(0, 10, 1))
  }
  else if (method == "rf"){
    tgrid <- expand.grid(.mtry = seq(10,100,10))
  }
  if (type == "class"){
    tcontrol <- trainControl(method="repeatedcv",
                             number= 4,
                             repeats = 4,
                             search="grid",
                             savePredictions ="all",
                             allowParallel = TRUE,
                             verboseIter=FALSE,
                             classProbs = TRUE,
                             sampling = "up")
    y <- y$class
  }
  else if (type == "regression"){
    tcontrol <- trainControl(method="repeatedcv",
                             number= 4,
                             repeats = 4,
                             search="grid",
                             savePredictions ="all",
                             allowParallel = TRUE,
                             verboseIter=FALSE)
    y <- y$aac
  }
  trainIndex <- createDataPartition(y, p = .8,
                                    list = TRUE,
                                    times = 25)
  pred_sample <- data.frame()
  per <- list()
  modRes <- list()
  output <- list()
  for(res in names(trainIndex))
  {
    trIndx <- trainIndex[[res]]
    tsIndx <- setdiff(1:nrow(x), trIndx)

    preProcValues <- preProcess(x[trIndx, ], method = c("center", "scale"))
    trainTransformed <- predict(preProcValues, x[trIndx, ])
    testTransformed <- predict(preProcValues, x[tsIndx, ])
    if (type == "class"){
      train_result_sample <- train(x=trainTransformed, y=y[trIndx],
                                   method=sprintf("%s", method),
                                   maximize = TRUE,
                                   tuneGrid=tgrid,
                                   trControl=tcontrol)
      train_result_sample$trainingData <- NULL
      predictedRes  <- predict(train_result_sample, testTransformed, type = c("raw"))
      predictedResProb  <- predict(train_result_sample, testTransformed, type = c("prob"))
      pred_sample_n <- data.frame(index = tsIndx,
                                  predict.class=predictedRes,
                                  predict.prob=predictedResProb,
                                  obs=y[tsIndx],
                                  resample=res)
      per <- confusionMatrix(data = as.factor(pred_sample_n$predict.class), reference = as.factor(pred_sample_n$obs))
      pred_sample <- rbind(pred_sample, pred_sample_n)
      metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))
      modRes[[res]] <- list("model"=train_result_sample, "preprocess" = preProcValues)
      output[[res]] <- list("prediction" = pred_sample, "stats" = per, "metadata" = metadata)
    } else if (type  == "regression"){
      if (ft > 0){
        featcor <- abs(apply(trainTransformed, 2, function(i) cor(i,y[trIndx])))
        featcor <- sort(featcor, decreasing = TRUE)
        featcor <- featcor[1:ft]
        trainTransformed <- trainTransformed[, names(featcor)]
        testTransformed <- testTransformed[, names(featcor)]
      } else if (var_count > 0){
        gene_vars <- abs(apply(x[trIndx,], 2, var))
        gene_vars <- sort(gene_vars, decreasing = TRUE)
        gene_vars <- gene_vars[1:var_count]
        trainTransformed <- trainTransformed[, names(gene_vars)]
        testTransformed <- testTransformed[, names(gene_vars)]
      }
      train_result_sample <- train(x=trainTransformed, y=y[trIndx],
                                   method=sprintf("%s", method),
                                   maximize = TRUE,
                                   tuneGrid=tgrid,
                                   trControl=tcontrol)
      train_result_sample$trainingData <- NULL
      pred_sample_n <- data.frame(index = tsIndx,
                                  pred=predict(train_result_sample, testTransformed),
                                  obs=y[tsIndx],
                                  resample=res)
      per <- postResample(pred = pred_sample_n$pred, obs = pred_sample_n$obs)
      pred_sample <- rbind(pred_sample, pred_sample_n)
      metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))
      modRes[[res]] <- list("model"=train_result_sample, "preprocess" = preProcValues)
      output[[res]] <- list("prediction" = pred_sample, "stats" = per, "metadata" = metadata)
    }
  }
  return(list("model" = modRes, "output" = output))
}
