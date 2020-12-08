library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(data.table)
source("./summarizeData.R")
source("./plot_data.R")
source("./computeInteractionMatrix.R")
GDSC2 <- readRDS("../data/GDSC2.rds")
drugname <- "Lapatinib"

generate_df <- function(pSet, mDataType, drug){
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drug, sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]
  
  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]
  
  #remove samples with no variance
  gene_vars <- apply(assay(df), 1, var)
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  gene_vars <- gene_vars[gene_vars>0]
  df <- df[names(gene_vars), ]
  
  return(df)
}

subset_by_feat <- function(df, drug, textmining, subset_size = 500){
  minedgenes <- readRDS(sprintf("./genes/%s.rds",toupper(drug)))
  #select mined genes
  if (textmining == TRUE){
    dfout = df[rowData(df)$gene_name %in% minedgenes$Symbol, ]
  } else {
    dfout= df[!(rowData(df)$gene_name %in% minedgenes$Symbol), ]
    dfout = dfout[1:subset_size]
  }
  
  #produce x and y, turn y to discrete
  x <- t(assay(dfout))
  y <- colData(dfout)$aac_recomputed
  cutoff <- callingWaterfall(y, "AUC")
  y2 <- ifelse(y >= cutoff,"sensitive","resistant")
  
  #print class counts for y
  counts <- table(y2)
  print(sprintf("Resistant: %d", counts[1]))
  print(sprintf("Sensitive: %d", counts[2]))
  
  output <- list("X"=x, "Y"=y2, "aac" = y)
  return(output)
}

#second function trains model based on x and y
trainmodel <- function(x,y,name,type){
  trainIndex <- createDataPartition(y, p = .8, 
                                    list = TRUE, 
                                    times = 10)
  
  tgrid <- expand.grid(alpha=seq(0, 1, 0.2),
                       lambda=seq(0, 10, 1))
  
  tcontrol <- trainControl(method="repeatedcv",
                           number= 4,
                           repeats = 4,
                           search="grid",
                           savePredictions ="all",
                           allowParallel = TRUE,
                           verboseIter=TRUE,
                           classProbs = TRUE)
  
  pred_sample <- data.frame()
  per <- list()
  
  for(res in names(trainIndex))
  {
    trIndx <- trainIndex[[res]]
    tsIndx <- setdiff(1:nrow(x), trIndx)
    
    #trainTransformed <- x[trIndx, ]
    #testTransformed <- x[tsIndx, ]
    preProcValues <- preProcess(x[trIndx, ], method = c("center", "scale"))
    trainTransformed <- predict(preProcValues, x[trIndx, ])
    testTransformed <- predict(preProcValues, x[tsIndx, ])
    
    if (type == "class"){
      train_result_sample <- train(x=trainTransformed, y=y[trIndx],
                                   method="glmnet",
                                   maximize = TRUE,
                                   tuneGrid=tgrid,
                                   trControl=tcontrol)
      train_result_sample$trainingData <- NULL
      pred_sample_n <- data.frame(index = tsIndx,
                                  predict=predict(train_result_sample, testTransformed),
                                  original=y[tsIndx], 
                                  resample=res)
      per[[res]] <- confusionMatrix(data = as.factor(pred_sample_n$predict), reference = as.factor(pred_sample_n$original))
      pred_sample <- rbind(pred_sample, pred_sample_n)
    }
    else if (type  == "regression"){
      train_result_sample <- train(x=trainTransformed, y=y[trIndx],
                                   method="glmnet",
                                   maximize = TRUE,
                                   tuneGrid=tgrid,
                                   trControl=tcontrol)
      train_result_sample$trainingData <- NULL
      pred_sample_n <- data.frame(index = tsIndx,
                                  predict=predict(train_result_sample, testTransformed),
                                  original=y[tsIndx], 
                                  resample=res)
      per[[res]] <- confusionMatrix(data = as.factor(pred_sample_n$predict), reference = as.factor(pred_sample_n$original))
      pred_sample <- rbind(pred_sample, pred_sample_n)
    }
  }
  
  #saveRDS(train_result_sample, sprintf("model_%s.rds", name))
  #saveRDS(pred_sample, sprintf("pred_result_%s.rds", name))
  metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))
  output <- list("prediction" = pred_sample, "stats" = per, "metadata" = metadata)
  return(output)
}

df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", drugname)

# for text mining genes
temp1 <- subset_by_feat(df, drugname, TRUE)
result1 <- trainmodel(temp1$X, temp1$Y, drugname, "class")
bwdata1 <- sapply(result1$stats, function(temp) temp$overall["Accuracy"])

# for top 500 genes
temp2 <- subset_by_feat(df, drugname, FALSE, 500)
result2 <- trainmodel(temp2$X, temp2$Y, drugname, "class")
bwdata2 <- sapply(result2$stats, function(temp) temp$overall["Accuracy"])

# for top 100 genes
temp3 <- subset_by_feat(df, drugname, FALSE, 100)
result3 <- trainmodel(temp3$X, temp3$Y, drugname, "class")
bwdata3 <- sapply(result3$stats, function(temp) temp$overall["Accuracy"])

# for top same number of genes
temp4 <- subset_by_feat(df, drugname, FALSE, result1$metadata$features)
result4 <- trainmodel(temp4$X, temp4$Y, drugname, "class")
bwdata4 <- sapply(result4$stats, function(temp) temp$overall["Accuracy"])

plotbw(drugname, bwdata1, bwdata2, bwdata3, bwdata4)
