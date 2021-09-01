library(tools)
library(Metrics)
library(SummarizedExperiment)
library(PharmacoGx)
library(caret)
library(psych)
library(stringi)
library(tidyverse)
library(tools)
library(data.table)
source("./train_functions.R")

pSet <- "CCLE"
testSet <- "GDSC"
problem <- "regression"
method <- "glmnet"
metric <- "COR"
genespath <- "./allgenes/"
type <- "train"

files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
#models <- list("tm", "500", "100", "ntm", "ft", "L1000","nL1000")
models <- list("tm", "tm100", "tm500", "tm1000")
if (pSet == "CCLE"){
  dataset <- readRDS("../data/CCLE_CTRPv2_solidTumor.rds")
  mDataType <- "rna"
  trainset <- "GDSC"
} else if (pSet == "GDSC") {
  dataset <- readRDS("../data/GDSC2.rds")
  mDataType <- "Kallisto_0.46.1.rnaseq"
  trainset <- "CCLE"
} else if (pSet == "gCSI") {
  dataset <- readRDS("../data/gCSI2.rds")
  mDataType <- "Kallisto_0.46.1.rnaseq" 
  trainset <- "CCLE"
}
ci <- cellInfo(dataset)
ci2 <- ci[!ci$tissueid %in% c("Lymphoid", "Myeloid"), ]
dataset <- subsetTo(dataset,cells = ci2$cellid)

jointdata <- list(tm = list(), L1000 = list(), ft = list(), ntm = list(), "100" = list(), "500" = list(), nL1000=list())
if (problem == "regression"){
  for (file in files){
    drugname <- file_path_sans_ext(file)
    print(drugname)
    if (type == "train"){
      tryCatch({
        for (model in models){
          print(model)
          data <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_%s.rds",pSet, problem, drugname, method,model))
          #cor <- sapply(data, function(temp) cor(temp$prediction$pred, temp$prediction$obs, method = "spearman"))
          cor <- sapply(data, function(temp) temp$stats["RMSE"])
          jointdata[[model]] <- cbind(jointdata[[model]], max(cor))
        }
      },
      error = function(e) { 
        print(e)
      })
    }
    else {
      tryCatch({
        for (model in models){
          print(model)
          data <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",testSet, problem, drugname, method,model))
          #cor <- sapply(data, function(temp) cor(temp$pred, temp$original, method = "spearman"))
          cor <- sapply(data, function(temp) temp$stats["RMSE"])
          jointdata[[model]] <- cbind(jointdata[[model]], max(cor))
        }
        }, 
        error = function(e) { 
          print(e)
          })
    }
  }
  saveRDS(jointdata, "../data/paired_cor.rds")
} else {
  for (file in files){
    drugname <- file_path_sans_ext(file)
    print(drugname)
    #df <- generate_df(dataset, mDataType, str_to_title(drugname))
    if (type == "train"){
      tryCatch({
        #temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
        for (model in models){
          print(model)
          data <- readRDS(sprintf("../train_output/%s/%s/output1/%s_%s_%s.rds",pSet, problem, drugname, method,model))
          if (metric == "AUC"){
            metricdata <- sapply(data, function(temp) auc(ifelse(temp$prediction$obs == "sensitive",1,0),temp$prediction$predict.prob.sensitive))
          } else if (metric == "Accuracy"){
            metricdata <- sapply(data, function(temp) temp$stats$overall["Accuracy"])
          } else if (metric == "BalancedAccuracy"){
            metricdata <- sapply(data, function(temp) temp$stats$byClass["Balanced Accuracy"])
          } else if (metric == "COR"){
            metricdata <- sapply(data, function(temp) cor(temp1$Y$aac[temp$prediction$index], temp$prediction$predict.prob.sensitive))
          }
          jointdata[[model]] <- cbind(jointdata[[model]], max(metricdata))
        }
      }, 
    error = function(e) { 
      print(e)
      })
    }
    else {
      tryCatch({
      for (model in models){
        print(model)
        #data <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",testSet, problem, drugname, method,model))
        data <- readRDS(sprintf("../test_output/%s/class/%s_%s_%s.rds", testSet, drugname, method, model))
        if (metric == "AUC"){
          metricdata <- sapply(data, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
        } else if (metric == "Accuracy"){
          metricdata <- sapply(data, function(temp) temp$stats$overall["Accuracy"])
        } else if (metric == "BalancedAccuracy"){
          metricdata <- sapply(data, function(temp) temp$stats$byClass["Balanced Accuracy"])
        } else if (metric == "COR"){
          metricdata <- sapply(data, function(temp) cor(temp$original, temp$pred$prob$sensitive))
        }
        jointdata[[model]] <- cbind(jointdata[[model]], max(metricdata))
      }
      }, 
      error = function(e) { 
        print(e)
      })
    }
  }
  saveRDS(jointdata, sprintf("../data/paired_%s.rds", metric))
}



