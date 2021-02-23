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
source("./mRMR.R")

set.seed(1)
GDSC2 <- readRDS("../data/GDSC2.rds")
ci <- cellInfo(GDSC2)
ci2 <- ci[!ci$tissueid %in% c("Lymphoid", "Myeloid"), ]
GDSC2 <- subsetTo(GDSC2,cells = ci2$cellid)

CCLE <- readRDS("../data/CCLE_CTRPv2_solidTumor.rds")
gene_count <- readRDS("../data/gene_count.rds")
L1000_gene_list <- read.delim("L1000_gene_list.txt")$Symbol


genespath <- "C:\\Users\\Grace Wu\\Documents\\text_features\\R\\genes\\"

drugname <- "Alisertib"
pSet <- "GDSC2"
method <- "glmnet"
problem <- "regression"

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

subset_by_feat <- function(df, drug, textmining = NULL , subset_size = 0, cutoff_method = "waterfall", L1000 = FALSE, ft = FALSE){
  drug <- "Alisertib"
  minedgenes <- readRDS(sprintf("./genes/%s.rds",toupper(drug)))
  #select mined genes
  if (L1000 == TRUE){
    dfout = df[rowData(df)$gene_name %in% L1000_gene_list, ]
  } else if (textmining == TRUE){
    dfout = df[rowData(df)$gene_name %in% minedgenes$Symbol, ]
    if (ft == TRUE && subset_size != 0){
      featcor <- apply(assay(dfout), 1, function(i) cor(i,dfout$aac_recomputed))
      featcor <- sort(featcor, decreasing = TRUE)
      featcor <- featcor[1:100]
      dfout <- dfout[names(featcor), ]
    }
  } else if (textmining == FALSE) {
    dfout= df[!(rowData(df)$gene_name %in% minedgenes$Symbol), ]
    if (ft == TRUE && subset_size != 0){
      featcor <- apply(assay(dfout), 1, function(i) cor(i,dfout$aac_recomputed))
      featcor <- sort(featcor, decreasing = TRUE)
      featcor <- featcor[1:100]
      dfout <- dfout[names(featcor), ]
    }
  }
  
  if (subset_size != 0){
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
trainmodel <- function(x,y,name,type, method, ft = -100, var = -100){
  trainIndex <- createDataPartition(y, p = .8,
                                    list = TRUE,
                                    times = 10)
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
  }
  else if (type == "regression"){
    tcontrol <- trainControl(method="repeatedcv",
                             number= 4,
                             repeats = 4,
                             search="grid",
                             savePredictions ="all",
                             allowParallel = TRUE,
                             verboseIter=FALSE)
  }
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
    set.seed(1)
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
      per <- confusionMatrix(data = as.factor(pred_sample_n$predict.class), reference = as.factor(pred_sample_n$original))
      pred_sample <- rbind(pred_sample, pred_sample_n)
      metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))##change this line
      modRes[[res]] <- list("model"=train_result_sample, "prediction" = pred_sample, 
                            "stats" = per, "metadata" = metadata)
    } else if (type  == "regression"){
      if (ft > 0){
        #trainTransformedtest = trainTransformed[,0:4000]
        #features <- .fs.mRMR(trainTransformedtest, y[trIndx], feature.count = 3897)
        #trainTransformed <- trainTransformed[, names(features)]
        #testTransformed <- testTransformed[, names(features)]
        featcor <- abs(apply(trainTransformed, 2, function(i) cor(i,y[trIndx])))
        featcor <- sort(featcor, decreasing = TRUE)
        featcor <- featcor[1:ft]
        trainTransformed <- trainTransformed[, names(featcor)]
        testTransformed <- testTransformed[, names(featcor)]
      } else if (var > 0){
        gene_vars <- abs(apply(trainTransformed, 2, function(i) var))
        gene_vars <- sort(gene_vars, decreasing = TRUE)
        gene_vars <- gene_vars[1:var]
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
      metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))##change this line
      modRes[[res]] <- list("model"=train_result_sample, "preprocess" = preProcValues)
      output[[res]] <- list("prediction" = pred_sample, "stats" = per, "metadata" = metadata)
    }
  }
  #saveRDS(train_result_sample, sprintf("model_%s.rds", name))
  #metadata <- list("samples" = nrow(x), "features" = ncol(x), "label" = table(y))
  return(list("model" = modRes, "output" = output))
}

files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)

for (file in files){
  drugname <- stri_sub(file, 1, -5)
  print(drugname)
  tryCatch({
    df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))

    print("text mining genes")
    temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
    for (i in colnames(temp1$X)){
      result1 <- trainmodel(temp1$X[,i], temp1$aac, drugname, problem ,method)
      saveRDS(result1$model, sprintf("../train_output/univariate/model/%s_%s_%s_%s_%s_tm.rds", pSet,drugname,method, problem, i))
      saveRDS(result1$output, sprintf("../train_output/univariate/output/%s_%s_%s_%s_%s_tm.rds", pSet,drugname,method, problem, i))
    }
    rm(temp1, result1)
    
  },error = function(e) {
    print("Drug not found in database.")})
}

drugname <- "Lapatinib"

univariate_train <- function(df, drugname, problem, method){
  temp4 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed", ft = TRUE, subset_size = 100)
  for (i in colnames(temp4$X)){
    result4 <- trainmodel(temp4$X, temp4$aac, drugname, problem, method)
    saveRDS(result4$model, sprintf("../train_output/univariate/model/Bortezomib_ntm/%s_%s_%s_%s_%s_ntm.rds", pSet,drugname,method, problem, i))
    saveRDS(result4$output, sprintf("../train_output/univariate/output/Bortezomib_ntm/%s_%s_%s_%s_%s_ntm.rds", pSet,drugname,method, problem, i))
  }
  rm(temp4, result4)
  temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed", ft = TRUE, subset_size = 100)
  for (i in colnames(temp1$X)){
    result1 <- trainmodel(temp1$X, temp1$aac, drugname, problem ,method)
    saveRDS(result1$model, sprintf("../train_output/univariate/model/Bortezomib_tm/%s_%s_%s_%s_%s_tm.rds", pSet,drugname,method, problem, i))
    saveRDS(result1$output, sprintf("../train_output/univariate/output/Bortezomib_tm/%s_%s_%s_%s_%s_tm.rds", pSet,drugname,method, problem, i))
  }
}

df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))
print("text mining genes")
temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
result1 <- trainmodel(temp1$X, temp1$aac, drugname, problem ,method)
#saveRDS(result1, sprintf("../train_output/%s_%s_%s_%s_tm.rds", pSet,drugname,method, problem))
saveRDS(result1, sprintf("../train_output/%s_%s_%s_%s_tm_ft.rds", pSet,drugname,method, problem))
rm(result1)

print("feature selection genes 100")
temp5 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
result5 <- trainmodel(temp5$X, temp5$aac, drugname, problem ,method, ft = 100)
#saveRDS(result5, sprintf("../train_output/%s_%s_%s_%s_ft.rds", pSet,drugname,method, problem))
saveRDS(result5, sprintf("../train_output/temp/%s_%s_%s_%s_ft_100.rds", pSet,drugname,method, problem))

print("feature selection genes 500")
temp5 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
result5 <- trainmodel(temp5$X, temp5$aac, drugname, problem ,method, ft = 500)
#saveRDS(result5, sprintf("../train_output/%s_%s_%s_%s_ft.rds", pSet,drugname,method, problem))
saveRDS(result5, sprintf("../train_output/temp/%s_%s_%s_%s_ft_500.rds", pSet,drugname,method, problem))

print("feature selection genes ft")
temp5 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
result5 <- trainmodel(temp5$X, temp5$aac, drugname, problem ,method)
#saveRDS(result5, sprintf("../train_output/%s_%s_%s_%s_ft.rds", pSet,drugname,method, problem))
saveRDS(result5, sprintf("../train_output/%s_%s_%s_%s_ft_100.rds", pSet,drugname,method, problem))

print("top 500 genes")
temp2 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
result2 <- trainmodel(temp2$X, temp2$aac, drugname, problem, method, var = 500)
saveRDS(result2, sprintf("../train_output/%s_%s_%s_%s_500.rds", pSet,drugname,method, problem))

print("top 100 genes")
temp3 <- subset_by_feat(df, drugname, FALSE,  cutoff_method = "fixed")
result3 <- trainmodel(temp3$X, temp3$aac, drugname, problem, method, var = 100)
saveRDS(result3, sprintf("../train_output/%s_%s_%s_%s_100.rds", pSet,drugname,method, problem))

tm <- readRDS(sprintf("../train_output/temp/%s_%s_%s_%s_tm_0.8.rds", pSet,drugname, method, problem))
print("not text mining")
temp4 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
result4 <- trainmodel(temp4$X, temp4$aac, drugname, problem, method, var = tm$Resample01$metadata$features)
saveRDS(result4, sprintf("../train_output/%s_%s_%s_%s_ntm.rds", pSet,drugname,method, problem))

print("L1000 genes")
temp6 <- subset_by_feat(df, drugname, cutoff_method = "fixed", L1000 = TRUE)
result6 <- trainmodel(temp6$X, temp6$aac, drugname, problem ,method)
saveRDS(result6, sprintf("../train_output/%s_%s_%s_%s_L1000.rds", pSet,drugname,method, problem))
