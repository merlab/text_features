#!/usr/bin/env Rscript
# source("./R/helper/summarizeData.R")
# source("./R/helper/computeInteractionMatrix.R")
# source("./R/helper/mRMR.R")
args <- commandArgs(trailingOnly = TRUE)
mlMethods <- c("rf")
geneSelectMethods <- c("text-mining")
# "var-500", "var-100", "L1000", "L1000-tm", "cor-500",
# # "glmnet",
library(future)
library(caret)
library(psych)
library(data.table)
library(stringi)
library(tools)
library(tidyverse)
#
library(doParallel)
library(foreach)
max.try <- as.numeric(args[3]) # 200

min.desirable.distance <- as.numeric(args[4])
max.desirable.distance <- as.numeric(args[5])

cl <- makePSOCKcluster(16)
registerDoParallel(cl)
#
calculateMetrics <- function(model, pred, obs) {
  model$perfMetric$pearsonCor <- cor(pred, obs, method = "pearson")
  model$perfMetric$spearmanCor <- cor(pred, obs, method = "spearman")
  model$perfMetric$kendallCor <- cor(pred, obs, method = "kendall")
  model$perfMetric$RMSE <- RMSE(pred = pred, obs = obs)
  model$perfMetric$MSE <- mean((obs - pred)^2)
  model$perfMetric$MAE <- MAE(pred = pred, obs = obs)
  return(model)
}

trainModel <- function(x, y, gdseX, gdseY,
                       mlMethod,
                       geneSelectMethod,
                       # tuningGrid,
                       corTargetList,
                       drug,
                       seed = 1,
                       min.desirable.distance = -0.01,
                       max.desirable.distance = 0.01) {
  # double check the harmonization
  sharedGenes <- intersect(colnames(x), colnames(gdseX))
  x <- x[, sharedGenes]
  gdseX <- gdseX[, sharedGenes]
  # set.seed(seed)
  # since this is regression we use this only
  tControl <- trainControl(
    allowParallel = TRUE,
    savePredictions = "none",
    returnResamp = "none",
    returnData = FALSE
    # NOTE:
  )

  indexes <- createDataPartition(
    seq_len(nrow(x)),
    p = 0.8, # runif(1, min = 0.4, max = 0.8), # .8,
    list = TRUE,
    times = 5
  )
  modelList <- list()

  bestModel <- list()
  for (i in names(indexes)[1]) {
    # modelList <- foreach(i == seq_len(indexes)) %dopar% {
    trIndx <- indexes[[i]]
    tsIndx <- setdiff(seq_len(nrow(x)), trIndx)
    subx <- x[trIndx, ]
    suby <- y[trIndx]

    # print(head(subx))
    gvar <- apply(subx, 2, var)
    gcor <- apply(subx, 2, function(a) {
      cor(a, suby, method = "pearson")
    })
    selectedGenes <- selectGenes(
      method = geneSelectMethod,
      genes = colnames(subx),
      geneCor = gcor,
      geneVar = gvar,
      drug = drug
    )
    selectedGenes <- intersect(selectedGenes, colnames(x))
    selectedGenes <- intersect(selectedGenes, colnames(gdseX))
    subx <- subx[, selectedGenes]
    gdseX <- gdseX[, selectedGenes]
    preProcValues <- preProcess(subx, method = c("center", "scale"))
    subXtrans <- predict(preProcValues, subx)
    model <- train(
      x = subXtrans,
      y = suby,
      method = mlMethod,
      metric = "RMSE",
      # tuningGrid = tuningGrid[(2 * j - 1):(2 * j), ],
      trControl = tControl
    )
    model$trainingData <- NULL
    model$finalModel$call <- NULL
    model$seed <- seed
    model$genes <- selectedGenes
    model$partition <- trIndx

    gdseXtrans <- predict(preProcValues, gdseX)
    pred <- predict(model, gdseXtrans)
    achievedCor <- cor(pred, gdseY, method = "pearson")
    targetCor <- corTargetList[[mlMethod]][drug, geneSelectMethod]
    model$distanceCor <- achievedCor - targetCor
    model <- calculateMetrics(model, pred = pred, obs = gdseY)
    model$idx <- indexes[[i]]
  }
  return(model)
}

selectGenes <- function(method, genes, geneVar, geneCor, drug) {
  l1000Genes <- read.table("./data/L1000_gene_list.txt", sep = "\t", header = TRUE, fill = TRUE)$Symbol
  l1000Genes <- na.omit(l1000Genes)
  txtMiningGenes <- readRDS(sprintf("./data/drug_text-features/%s.rds", drug))$Symbol
  geneVar <- na.omit(geneVar)
  geneCor <- na.omit(geneCor)
  sortedGeneVar <- sort(geneVar, decreasing = TRUE)
  sortedGeneCor <- sort(abs(geneCor), decreasing = TRUE)
  if (method == "var-100") {
    sortedGenes <- names(sortedGeneVar)
    genes <- intersect(genes, sortedGenes[seq_len(100)])
  }
  if (method == "var-500") {
    sortedGenes <- names(sortedGeneVar)
    genes <- intersect(genes, sortedGenes[seq_len(500)])
  }
  if (method == "L1000-tm") {
    # genes <- intersect(genes, txtMiningGenes)
    genes <- intersect(genes, l1000Genes)
    genes <- genes[!genes %in% txtMiningGenes]
  }
  if (method == "L1000") {
    genes <- intersect(genes, l1000Genes)
  }
  if (method == "cor-500") {
    sortedGenes <- names(sortedGeneCor)
    genes <- intersect(genes, sortedGenes[seq_len(500)])
  }
  if (method == "cor-100") {
    sortedGenes <- names(sortedGeneCor)
    genes <- intersect(genes, sortedGenes[seq_len(100)])
  }
  if (method == "text-mining") {
    genes <- intersect(genes, txtMiningGenes)
  }
  return(genes)
}
procRNA <- function(mat) {
  rvar <- apply(mat, 1, var)
  mat <- mat[!is.na(rvar), ]
  cvar <- apply(mat, 2, var)
  mat <- mat[, !is.na(cvar)]
  return(mat)
}

trainIn <- readRDS("./result/ccle_ctrpv2_rnaseq.rds")
trainOut <- readRDS("./result/ccle_ctrpv2_aac.rds")
targetIn <- readRDS("./result/gdse_rnaseq.rds")
targetOut <- readRDS("./result/gdse_aac.rds")

tGrid <- list()

# harmonize columns
drugs <- intersect(colnames(trainOut), colnames(targetOut))
genes <- intersect(colnames(trainIn), colnames(targetIn))
trainIn <- trainIn[, genes]
targetIn <- targetIn[, genes]

trainOut <- trainOut[, drugs]
targetOut <- targetOut[, drugs]
#
trainIn <- procRNA(trainIn)
targetIn <- procRNA(targetIn)
#
genes <- intersect(colnames(trainIn), colnames(targetIn))
trainIn <- trainIn[, genes]
targetIn <- targetIn[, genes]
#
# harmonize rows
trainSamples <- intersect(rownames(trainIn), rownames(trainOut))
targetSamples <- intersect(rownames(targetIn), rownames(targetOut))
trainIn <- trainIn[trainSamples, ]
trainOut <- trainOut[trainSamples, ]
#
targetIn <- targetIn[targetSamples, ]
targetOut <- targetOut[targetSamples, ]
#

corTargetList <- list(
  "glmnet" = readRDS("./result/train_test/test_ElasticNet_Reg.rds"),
  "rf" = readRDS("./result/train_test/test_RandomForest_Reg.rds")
)



for (mlMethod in mlMethods) {
  for (geneSelectMethod in geneSelectMethods) {
    print(paste(mlMethod, geneSelectMethod))
    time1 <- Sys.time()
    for (drug in rev(drugs)) {
      name <- paste0(paste(drug, mlMethod, geneSelectMethod, sep = "-"), ".rds")
      trainX <- trainIn
      cvar <- apply(trainX, 2, var)
      targetX <- targetIn
      #
      sharedGenes <- intersect(colnames(trainX), colnames(targetX))
      targetX <- targetX[, sharedGenes]
      trainX <- trainX[, sharedGenes]
      #
      trainY <- as.vector(trainOut[, drug])
      trainY <- setNames(trainY, rownames(trainOut))
      targetY <- as.vector(targetOut[, drug])
      targetY <- setNames(targetY, rownames(targetOut))
      #
      trainY <- trainY[!is.na(trainY)]
      targetY <- targetY[!is.na(targetY)]
      trainX <- trainX[names(trainY), ]
      targetX <- targetX[names(targetY), ]
      #
      if (length(trainY) == 0 || length(targetY) == 0 ||
        nrow(trainX) == 0 || ncol(trainX) == 0 ||
        nrow(targetX) == 0 || ncol(targetX) == 0 ||
        var(trainY) == 0
      ) {
        next()
      }
      trainCut <- 0.75
      set.seed(1)
      rows <- seq_len(nrow(trainX))
      idx <- sample(rows, floor(nrow(trainX) * trainCut), replace = FALSE, prob = NULL) # runif(0.3, 0.9)
      model <- trainModel(
        x = trainX[idx, ],
        y = trainY[idx],
        gdseX = targetX,
        gdseY = targetY,
        drug = drug,
        mlMethod = mlMethod,
        geneSelectMethod = geneSelectMethod,
        corTargetList = corTargetList,
        seed = 1,
        min.desirable.distance = min.desirable.distance,
        max.desirable.distance = max.desirable.distance
      )
      print(Sys.time() - time1)
    }
    print(Sys.time() - time1)
  }
}
stopCluster(cl)
print("done")
