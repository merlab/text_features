drugname <- "Erlotinib"
pSet <- "GDSC2"
metric <- "RMSE"
method <- "glmnet"
problem <- "regression"

CCLE <- readRDS("../data/CCLE_CTRPv2_solidTumor.rds")
ccle_df <- generate_df(CCLE, "rna", str_to_title(drugname))

predictmodel <- function(pSet, drugname, method, problem){
  tm <- readRDS(sprintf("../train_output/%s_%s_%s_%s_tm.rds", pSet,drugname, method, problem))
  return(0)
  top500 <- readRDS(sprintf("../train_output/%s_%s_%s_%s_500.rds", pSet,drugname, method, problem))
  top100 <- readRDS(sprintf("../train_output/%s_%s_%s_%s_100.rds", pSet,drugname, method, problem))
  ntm <- readRDS(sprintf("../train_output/%s_%s_%s_%s_ntm.rds", pSet,drugname, method, problem))
  ft <- readRDS(sprintf("../train_output/%s_%s_%s_%s_ft.rds", pSet,drugname, method, problem))
  L1000 <- readRDS(sprintf("../train_output/%s_%s_%s_%s_L1000.rds", pSet,drugname, method, problem))
  models <- list(tm, top500, top100, ntm, ft, L1000)
  return(0)
  for (model in models){
    fe <- sapply(model, function(temp) temp$model$finalModel$xNames)
    valid <- 1
    for (i in names(fe[1,])){
      if (FALSE == all (fe[,i] %in% rownames(ccle_df))){
        print(sprintf("%s",i))
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
      print(sprintf("%s",i))
      ccle_test <- t(assay(ccle_df)[fe[,i],])
      temppredict <- predict(model[[i]]["model"], newdata = ccle_test)
      pred <- temppredict
      stats <- postResample(pred = as.numeric(unlist(temppredict)), obs = ccle_df$aac_recomputed)
      modRes[[i]] <- list("pred" = pred, "stats" = stats)
    }
    saveRDS(modRes, sprintf("../test_output/%s_%s_%s_%s_%s.rds", pSet,drugname,method, problem,model))
  }
  return(1)
}

temp <- predictmodel(GDSC2, drugname, method, problem)

