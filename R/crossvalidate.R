drugname <- "Dasatinib"
pSet <- "GDSC2"
metric <- "RMSE"
method <- "glmnet"
problem <- "regression"

generate_testing_df <- function(pSet, mDataType, drug){
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drug, sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]
  
  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]
  
  return(df)
}

CCLE <- readRDS("../data/CCLE_CTRPv2_solidTumor.rds")

predictmodel <- function(pSet, drugname, method, problem){
  models <- list("tm", "500", "100", "ntm", "ft_tm", "L1000")
  for (model in models){
    print(model)
    data <- readRDS(sprintf("../train_output/%s_%s_%s_%s_%s.rds", pSet,drugname, method, problem,model))
    fe <- sapply(data, function(temp) temp$model$finalModel$xNames)
    valid <- 1
    for (i in names(fe[1,])){
      if (FALSE == all (fe[,i] %in% rownames(ccle_df))){
        print(fe[which(!fe[,i] %in% rownames(ccle_df))])
        valid <- 0
      }
    }
    if (valid == 0){
      print("Features do not match - Cannot predict")
      return(0)
    }
    modRes <- list()
    for (i in names(fe[1,])){
      ccle_test <- t(assay(ccle_df)[fe[,i],])
      temppredict <- predict(data[[i]]["model"], newdata = ccle_test)
      pred <- temppredict
      stats <- postResample(pred = as.numeric(unlist(temppredict)), obs = ccle_df$aac_recomputed)
      modRes[[i]] <- list("pred" = pred, "stats" = stats)
    }
    saveRDS(modRes, sprintf("../test_output/%s_%s_%s_%s_%s.rds", pSet,drugname,method, problem,model))
  }
  return(1)
}

ccle_df <- generate_testing_df(CCLE, "rna", str_to_title(drugname))
temp <- predictmodel(pSet, drugname, method, problem)

files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)

for (file in files){
  drugname <- stri_sub(file, 1, -5)
  print(drugname)
  tryCatch({
    ccle_df <- generate_testing_df(CCLE, "rna", str_to_title(drugname))
    print("predicting")
    temp <- predictmodel(pSet, drugname, method, problem)
  },error = function(e) {
    print("Drug not found in database.")})
}

