library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(data.table)
library(stringi)
library(tools)
library(tidyverse)
library(doParallel)
source("./train_functions.R")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
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

set.seed(1)

cl <- makePSOCKcluster(36)
registerDoParallel(cl)

# Read Datasets
L1000_gene_list <- read.delim("L1000_gene_list.txt")$Symbol

# Parameters 
genespath <- "./genes/"

if (pSet == "CCLE"){
  dataset <- readRDS("../data/CCLE_CTRPv2_solidTumor.rds")
  mDataType <- "rna"
} else if (pSet == "GDSC") {
  GDSC2 <- readRDS("../data/GDSC2.rds")
  ci <- cellInfo(GDSC2)
  ci2 <- ci[!ci$tissueid %in% c("Lymphoid", "Myeloid"), ]
  GDSC2 <- subsetTo(GDSC2,cells = ci2$cellid)
  dataset <- GDSC2
  mDataType <- "Kallisto_0.46.1.rnaseq"
}

if (is.null(drugname)) {
  files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
} else {
  files = list(drugname)
}

for (file in files){
  drugname <- file_path_sans_ext(file)
  print(drugname)
  df <- generate_df(dataset, mDataType, str_to_title(drugname))
  
  print("text mining genes")
  temp1 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
  result1 <- trainmodel(temp1$X, temp1$Y, drugname, problem ,method)
  saveRDS(result1$model, sprintf("../train_output/%s/%s/model/%s_%s_tm.rds", pSet,problem, drugname,method))
  saveRDS(result1$output, sprintf("../train_output/%s/%s/output/%s_%s_tm.rds", pSet,problem,drugname,method))
  tm <- ncol(temp1$X)
  rm(result1, temp1)
  
  print("not text mining")
  temp4 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
  result4 <- trainmodel(temp4$X, temp4$Y, drugname, problem, method, var_count = tm)
  saveRDS(result4$model, sprintf("../train_output/%s/%s/model/%s_%s_ntm.rds", pSet,problem,drugname,method))
  saveRDS(result4$output, sprintf("../train_output/%s/%s/output/%s_%s_ntm.rds", pSet,problem,drugname,method))
  rm(result4, temp4)
  
  print("feature selection genes ft")
  temp5 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
  result5 <- trainmodel(temp5$X, temp5$Y, drugname, problem ,method, ft = tm)
  saveRDS(result5$model, sprintf("../train_output/%s/%s/model/%s_%s_ft.rds", pSet,problem,drugname,method))
  saveRDS(result5$output, sprintf("../train_output/%s/%s/output/%s_%s_ft.rds", pSet,problem,drugname,method))
  rm(result5, temp5)
  
  print("top 500 genes")
  temp2 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
  result2 <- trainmodel(temp2$X, temp2$Y, drugname, problem, method, var_count = 500)
  saveRDS(result2$model, sprintf("../train_output/%s/%s/model/%s_%s_500.rds", pSet,problem,drugname,method))
  saveRDS(result2$output, sprintf("../train_output/%s/%s/output/%s_%s_500.rds", pSet,problem,drugname,method))
  rm(result2, temp2)
  
  print("top 100 genes")
  temp3 <- subset_by_feat(df, drugname, FALSE,  cutoff_method = "fixed")
  result3 <- trainmodel(temp3$X, temp3$Y, drugname, problem, method, var_count = 100)
  saveRDS(result3$model, sprintf("../train_output/%s/%s/model/%s_%s_100.rds", pSet,problem,drugname,method))
  saveRDS(result3$output, sprintf("../train_output/%s/%s/output/%s_%s_100.rds", pSet,problem,drugname,method))
  rm(result3, temp3)
  
  
  print("L1000 genes")
  temp6 <- subset_by_feat(df, drugname, cutoff_method = "fixed", L1000 = TRUE)
  result6 <- trainmodel(temp6$X, temp6$Y, drugname, problem ,method)
  saveRDS(result6$model, sprintf("../train_output/%s/%s/model/%s_%s_L1000.rds", pSet,problem,drugname,method))
  saveRDS(result6$output, sprintf("../train_output/%s/%s/output/%s_%s_L1000.rds", pSet,problem,drugname,method))
  rm(result6, temp6)
  rm(df)
}

stopCluster(cl)