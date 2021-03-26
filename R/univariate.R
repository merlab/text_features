library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(data.table)
library(stringi)
library(tools)
library(tidyverse)
source("./train_functions.R")

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 1) {
  stop("Please supply arguments: pSet, drugname (optional)", call.=FALSE)
} else if (length(args)==2) {
  pSet <- args[1]
  drugname <- args[2]
} else if (length(args)==1) {
  pSet <- args[1]
  drugname <- NULL
}

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

genespath <- "./genes/"
if (is.null(drugname)) {
  files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
} else {
  files = list(drugname)
}

for (file in files){
  drugname <- file_path_sans_ext(file)
  print(drugname)
  tryCatch({
    df <- generate_df(dataset, mDataType, str_to_title(drugname))
    
    temp2 <- subset_by_feat(df, drugname, FALSE, cutoff_method = "fixed")
    featcor <- apply(temp2$X, 2, function(i) cor(i,temp2$Y$aac))
    saveRDS(featcor, sprintf("../train_output/univariate/%s/%s_ntm.rds", pSet,drugname))
    rm(featcor, temp2)
    
    temp2 <- subset_by_feat(df, drugname, TRUE, cutoff_method = "fixed")
    featcor <- apply(temp2$X, 2, function(i) cor(i,temp2$Y$aac))
    saveRDS(featcor, sprintf("../train_output/univariate/%s/%s_tm.rds", pSet,drugname))
    rm(featcor, temp2, df)
  },error = function(e) {
    print("Drug not found in database.")})
}