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

drugname <- "Alisertib"

for (file in files){
  drugname <- file_path_sans_ext(file)
  print(drugname)
  tryCatch({
    
    df <- generate_df(CCLE, mDataType, str_to_title(drugname))
    
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

sample_count_ccle <- data.frame(matrix(ncol = 2, nrow = 0))
sample_count_gdsc <- data.frame(matrix(ncol = 2, nrow = 0))

for (drug in gene_count$name){
  print(drug)
  tryCatch({
    df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))
    sample_count_gdsc <- rbind(sample_count_gdsc, c(drugname, length(df$aac_recomputed)))
    rm(df)
    df <- generate_df(CCLE, "rna", str_to_title(drugname))
    sample_count_ccle <- rbind(sample_count_ccle, c(drugname, length(df$aac_recomputed)))
    rm(df)
  },error = function(e) {
    print("Drug not found in database.")})
}
names <- c("name", "count")
colnames(sample_count_ccle) <- names
colnames(sample_count_gdsc) <- names
saveRDS(sample_count_ccle, paste("../data/", "sample_count_ccle.rds"))
saveRDS(sample_count_gdsc, paste("../data/", "sample_count_gdsc.rds"))
