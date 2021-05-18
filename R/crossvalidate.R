library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(data.table)
library(stringi)
library(tools)
library(tidyverse)
source("./summarizeData.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Please supply arguments: pSet, method, problem, drugname (optional)", call.=FALSE)
} else if (length(args)==4) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  drugname <- args[4]
} else if (length(args)==3) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  drugname <- NULL
}

print(sprintf("pSet: %s, method: %s, problem: %s", pSet, method, problem))

genespath <- "./genes/"

if (pSet == "CCLE"){
  dataset <- readRDS("../data/CCLE_CTRPv2_solidTumor.rds")
  mDataType <- "rna"
  trainset <- "GDSC"
  
} else if (pSet == "GDSC") {
  GDSC2 <- readRDS("../data/GDSC2.rds")
  ci <- cellInfo(GDSC2)
  ci2 <- ci[!ci$tissueid %in% c("Lymphoid", "Myeloid"), ]
  GDSC2 <- subsetTo(GDSC2,cells = ci2$cellid)
  dataset <- GDSC2
  mDataType <- "Kallisto_0.46.1.rnaseq"
  trainset <- "CCLE"
}

generate_testing_df <- function(pSet, mDataType, drug){
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drug, sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]
  
  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]
  
  return(df)
}

predictmodel <- function(pSet, trainset, drugname, method, problem){
  models <- list("tm", "500", "100", "ntm", "ft", "L1000")
  for (model in models){
    print(model)
    data <- readRDS(sprintf("../train_output/%s/%s/model/%s_%s_%s.rds",trainset, problem, drugname, method,model))
    fe <- sapply(data, function(temp) temp$model$finalModel$xNames)
    valid <- 1
    for (i in names(fe[1,])){
      if (FALSE == all (fe[,i] %in% rownames(df))){
        print(fe[which(!fe[,i] %in% rownames(df))])
        valid <- 0
      }
    }
    if (valid == 0){
      print("Features do not match - Cannot predict")
      return(0)
    }
    modRes <- list()
    for (i in names(fe[1,])){
      test<- t(assay(df)[fe[,i],])
      #preProcValues <- preProcess(test, method = c("center", "scale"))
      #test <- predict(preProcValues, test)
      test <- predict(data[[i]]$preprocess, test)
      if (problem == "regression"){
        temppredict <- predict(data[[i]]["model"], newdata = test)
        pred <- temppredict$model
        stats <- postResample(pred = as.numeric(unlist(temppredict)), obs = df$aac_recomputed)
      } else if (problem == "class") {
        y2 <- ifelse(df$aac_recomputed >= 0.2,"sensitive","resistant")
        probpredict <- predict(data[[i]]["model"], newdata = test, type = c("prob"))
        rawpredict <- predict(data[[i]]["model"], newdata = test)
        pred <- list("class" = rawpredict$model, "prob" = probpredict$model, "obs" = as.vector(y2))
        stats <- confusionMatrix(data = as.factor(unlist(rawpredict)), reference = as.factor(y2))
      }
      modRes[[i]] <- list("pred" = pred, "original" = df$aac_recomputed, "stats" = stats)
    }
    saveRDS(modRes, sprintf("../test_output/%s/%s/%s_%s_%s.rds", pSet, problem, drugname,method,model))
  }
  return(1)
}

if (is.null(drugname)){
  files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
  
  for (file in files){
    drugname <- stri_sub(file, 1, -5)
    print(drugname)
    df <- generate_testing_df(dataset, mDataType, str_to_title(drugname))
    print("predicting")
    temp <- predictmodel(pSet, trainset, drugname, method, problem)
  }
} else {
  df <- generate_testing_df(dataset, mDataType, str_to_title(drugname))
  predictmodel(pSet, trainset, drugname, method, problem)
}

