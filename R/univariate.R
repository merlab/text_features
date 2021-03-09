library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(data.table)
library(stringi)
library(tools)
library(tidyverse)
source("./train_test.R")
source("./train_test.R")

files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)

drugname <- "Erlotinib"

for (file in files){
  drugname <- file_path_sans_ext(file)
  print(drugname)
  tryCatch({
    
    df <- generate_df(CCLE, mDataType, str_to_title(drugname))
    
    temp2 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
    featcor <- apply(temp2$X, 2, function(i) cor(i,temp2$Y$aac))
    saveRDS(featcor, sprintf("../train_output/univariate/%s/%s_ntm.rds", pSet,drugname))
    
    temp2 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
    featcor <- apply(temp2$X, 2, function(i) cor(i,temp2$Y$aac))
    saveRDS(featcor, sprintf("../train_output/univariate/%s/%s_tm.rds", pSet,drugname))
    
  },error = function(e) {
    print("Drug not found in database.")})
}
