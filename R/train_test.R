library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(data.table)
library(stringi)
source("./summarizeData.R")
source("./computeInteractionMatrix.R")

GDSC2 <- readRDS("../data/GDSC2.rds")
ci <- cellInfo(GDSC2)
ci2 <- ci[!ci$tissueid %in% c("Lymphoid", "Myeloid"), ]
GDSC2 <- subsetTo(GDSC2,cells = ci2$cellid)

genespath <- "C:\\Users\\Grace Wu\\text_features\\R\\genes\\"

drugname <- "Dasatinib"
pSet <- "GDSC2"
method <- "rf"

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

  ##------code to remove features that have low variance ------
  ##atleast 15% should be non zero
  cutoff <- floor(ncol(df)*0.15)
  unique_count <- apply(assay(df), 1, function(v) length(unique(round(v, 4))))
  gn2take <- unique_count > cutoff
  df <- df[gn2take, ]
  ##------------------------------------------------------
  return(df)
}

subset_by_feat <- function(df, drug, textmining, subset_size = 500, cutoff_method = "waterfall"){
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

  output <- list("X"=x, "Y"=y2, "aac" = y)
  return(output)
}

#second function trains model based on x and y
trainmodel <- function(x,y,name,type, method){
  trainIndex <- createDataPartition(y, p = .8,
                                    list = TRUE,
                                    times = 10)

  if(method == "glmnet"){
    tgrid <- expand.grid(alpha=seq(0, 1, 0.2),
                         lambda=seq(0, 10, 1))
  }
  else if (method == "rf"){
    tgrid <- expand.grid(.mtry = (1:5))
  }

  tcontrol <- trainControl(method="repeatedcv",
                           number= 4,
                           repeats = 4,
                           search="grid",
                           savePredictions ="all",
                           allowParallel = TRUE,
                           verboseIter=FALSE,
                           classProbs = TRUE,
                           sampling = "up")

  pred_sample <- data.frame()
  per <- list()

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
                                   ntree = 100,
                                   trControl=tcontrol)
      train_result_sample$trainingData <- NULL
      predictedRes  <- predict(train_result_sample, testTransformed, type = c("raw"))
      predictedResProb  <- predict(train_result_sample, testTransformed, type = c("prob"))
      pred_sample_n <- data.frame(index = tsIndx,
                                  predict.class=predictedRes,
                                  predict.prob=predictedResProb,
                                  original=y[tsIndx],
                                  resample=res)
      per[[res]] <- confusionMatrix(data = as.factor(pred_sample_n$predict.class), reference = as.factor(pred_sample_n$original))
      pred_sample <- rbind(pred_sample, pred_sample_n)
    } else if (type  == "regression"){
      train_result_sample <- train(x=trainTransformed, y=y[trIndx],
                                   method=sprintf("%s", method),
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
  metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))
  output <- list("prediction" = pred_sample, "stats" = per, "metadata" = metadata)
  return(output)
}

#ifelse(!dir.exists(file.path("../train_output/", sprintf("%s",drugname))), dir.create(file.path("../train_output/", sprintf("%s",drugname))), FALSE)

files <- list.files(path=genespath, pattern="*.rds", full.names=FALSE, recursive=FALSE)

for (file in files){
  drugname <- stri_sub(file, 2, -5)
  print(drugname)
  tryCatch({
    df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))

    print("text mining genes")
    temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
    result1 <- trainmodel(temp1$X, temp1$Y, drugname, "class",method)
    saveRDS(result1, sprintf("../train_output/%s_%s_%s_tm.rds", pSet,drugname,method))

    print("top 500 genes")
    temp2 <- subset_by_feat(df, drugname, FALSE, 500,  cutoff_method = "fixed")
    result2 <- trainmodel(temp2$X, temp2$Y, drugname, "class", method)
    saveRDS(result2, sprintf("../train_output/%s_%s_%s_500.rds", pSet,drugname,method))

    print("top 100 genes")
    temp3 <- subset_by_feat(df, drugname, FALSE, 100,  cutoff_method = "fixed")
    result3 <- trainmodel(temp3$X, temp3$Y, drugname, "class", method)
    saveRDS(result3, sprintf("../train_output/%s_%s_%s_100.rds", pSet,drugname,method))

    print("not text mining genes")
    temp4 <- subset_by_feat(df, drugname, FALSE, result1$metadata$features, cutoff_method = "fixed")
    result4 <- trainmodel(temp4$X, temp4$Y, drugname, "class", method)
    saveRDS(result4, sprintf("../train_output/%s_%s_%s_ntm.rds", pSet,drugname,method))
  },error = function(e) {
    print("Drug not found in database.")})
}


df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))

print("text mining genes")
temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
result1 <- trainmodel(temp1$X, temp1$Y, drugname, "class",method)
saveRDS(result1, sprintf("../train_output/%s_%s_%s_tm.rds", pSet,drugname,method))

print("top 500 genes")
temp2 <- subset_by_feat(df, drugname, FALSE, 500,  cutoff_method = "fixed")
result2 <- trainmodel(temp2$X, temp2$Y, drugname, "class", method)
saveRDS(result2, sprintf("../train_output/%s_%s_%s_500.rds", pSet,drugname,method))

print("top 100 genes")
temp3 <- subset_by_feat(df, drugname, FALSE, 100,  cutoff_method = "fixed")
result3 <- trainmodel(temp3$X, temp3$Y, drugname, "class", method)
saveRDS(result3, sprintf("../train_output/%s_%s_%s_100.rds", pSet,drugname,method))

print("not text mining")
temp4 <- subset_by_feat(df, drugname, FALSE, result1$metadata$features, cutoff_method = "fixed")
result4 <- trainmodel(temp4$X, temp4$Y, drugname, "class", method)
saveRDS(result4, sprintf("../train_output/%s_%s_%s_ntm.rds", pSet,drugname,method))
