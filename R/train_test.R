library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
source("./summarizeData.R")
GDSC2 <- readRDS("../data/GDSC2.rds")

#first function takes pset, mDataType and drug name to return x and y variables for training
dataset_preprocess <- function(pSet, mDataType, drug){
  
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drug, sensitivity.measure="aac_recomputed")
  df = summarizeData(pSet=GDSC2, mDataType = "Kallisto_0.46.1.rnaseq", drug = "Erlotinib", sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]
  
  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]
  
  #select mined genes
  minedgenes <- readRDS("./ERLOTINIB.rds")
  df = df[rowData(df)$gene_name <= minedgenes$Symbol, ]

  #remove samples with no variance
  gene_vars <- apply(assay(df), 1, var)
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  gene_vars <- gene_vars[gene_vars>0]
  df <- df[names(gene_vars), ]
  df = df[1:500]
  
  #produce x and y, turn y to discrete
  x <- t(assay(df))
  y <- colData(df)$aac_recomputed
  y <- ifelse(y >= 0.1,"sensitive","resistant")
  
  #print class counts for y
  counts <- table(y)
  sprintf("Resistant: %d", counts[1])
  sprintf("Sensitive: %d", counts[2])
  
  output <- list(x,y)
  return(output)
}

#second function trains model based on x and y
trainmodel <- function(x,y,name){
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
  
  saveRDS(train_result_sample, sprintf("model_%s.rds", name))
  saveRDS(pred_sample, sprintf("pred_result_%s.rds", name))
  metadata <- list(table(y), dim(x))
  output <- list(train_result_sample, pred_sample, metadata)
  return(output)
}

temp <- dataset_preprocess(GDSC2, "Kallisto_0.46.1.rnaseq", "Erlotinib")
result <- trainmodel(output[1][[1]], output[2][[1]], "erlotinib")
bwplot(result [[1]]$resample$Accuracy,  xlab="Accuracy", ylab = "Erlotinib", main = "Text Mining - Accuracy by Drug")


