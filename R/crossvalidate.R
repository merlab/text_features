args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("Please supply arguments: pSet, method, problem, drugname (optional)", call.=FALSE)
} else if (length(args)==4) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  drugname <- args[4]
  print("here1")
} else if (length(args)==3) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  drugname <- NULL
  print("here2")
}

#drugname <- "Axitinib"
#pSet <- "GDSC2"
#method <- "glmnet"
#problem <- "regression"
genespath <- "C:\\Users\\Grace Wu\\Documents\\text_features\\R\\genes\\"

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


generate_testing_df <- function(pSet, mDataType, drug){
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drug, sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]
  
  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]
  
  return(df)
}

predictmodel <- function(pSet, drugname, method, problem){
  models <- list("tm", "500", "100", "ntm", "ft", "L1000")
  for (model in models){
    print(model)
    data <- readRDS(sprintf("../train_output/CCLE/regression/model/%s_%s_%s_%s.rds",drugname, method, problem,model))
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
      test <- t(assay(df)[fe[,i],])
      temppredict <- predict(data[[i]]["model"], newdata = test)
      pred <- temppredict$model
      stats <- postResample(pred = as.numeric(unlist(temppredict)), obs = df$aac_recomputed)
      modRes[[i]] <- list("pred" = pred, "original" = df$aac_recomputed, "stats" = stats)
    }
    saveRDS(modRes, sprintf("../test_output/GDSC2/%s_%s_%s_%s.rds", drugname,method, problem,model))
  }
  return(1)
}

drugname <- "Alpelisib"

if (is.null(drugname)){
  files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
  
  for (file in files){
    drugname <- stri_sub(file, 1, -5)
    print(drugname)
    tryCatch({
      df <- generate_testing_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))
      print("predicting")
      temp <- predictmodel(pSet, drugname, method, problem)
    },error = function(e) {
      print("Drug not found in database.")})
  }
} else {
  df <- generate_testing_df(GDSC2, "Kallisto_0.46.1.rnaseq", str_to_title(drugname))
  temp <- predictmodel(pSet, drugname, method, problem)
}
