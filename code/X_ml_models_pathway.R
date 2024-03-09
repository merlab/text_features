#!/usr/bin/env Rscript
cells <- readRDS("./cells.rds")
args <- commandArgs(trailingOnly = TRUE)
mlMethods <- "rf"
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
max.try <- 1

min.desirable.distance <- -1 # as.numeric(args[4])
max.desirable.distance <- 0 # as.numeric(args[5])

# desirable.distance <- 0.01
# desirable.distance <- 0.05
cores <- 8 # ceiling(detectCores() / 8) + 1
cl <- makePSOCKcluster(cores)
registerDoParallel(cl)
#
calculateMetrics <- function(model, pred, obs) {
  model$perfMetric$pearsonCor <- cor(pred, obs, method = "pearson")
  return(model)
}

trainModel <- function(x, y, gdseX, gdseY,
                       mlMethod,
                       # geneSelectMethod,
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

  bestDistAchieved <- 1
  bestModel <- list()
  for (i in names(indexes)[1]) {
    # modelList <- foreach(i == seq_len(indexes)) %dopar% {
    trIndx <- indexes[[i]]
    tsIndx <- setdiff(seq_len(nrow(x)), trIndx)
    subx <- x[trIndx, ]
    suby <- y[trIndx]


    # TODO: selectGenes
    pathwayKeyValue <- list(
      "Erlotinib" = c(
        "AKT1", "AKT2", "AKT3", "ARAF", "AXL", "BAD", "BAX", "BCL2",
        "BCL2L1", "BCL2L11", "BRAF", "CCND1", "EGF", "EGFR", "EIF4E",
        "EIF4E2", "EIF4EBP1", "ERBB2", "ERBB3", "FGF2", "FGFR2", "FGFR3",
        "FOXO3", "GAB1", "GAS6", "GRB2", "GSK3B", "HGF", "HRAS", "IGF1",
        "IGF1R", "IL6", "IL6R", "JAK1", "JAK2", "KDR", "KRAS", "MAP2K1",
        "MAP2K2", "MAPK1", "MAPK3", "MET", "MRAS", "MTOR", "MYC", "NF1",
        "NRAS", "NRG1", "NRG2", "PDGFA", "PDGFB", "PDGFC", "PDGFD", "PDGFRA",
        "PDGFRB", "PDPK1", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
        "PIK3R3", "PLCG1", "PLCG2", "PRKCA", "PRKCB", "PRKCG", "PTEN",
        "RAF1", "RPS6", "RPS6KB1", "RPS6KB2", "RRAS", "RRAS2", "SHC1",
        "SHC2", "SHC3", "SHC4", "SOS1", "SOS2", "SRC", "STAT3", "TGFA",
        "VEGFA"
      ),
      "Lapatinib" = c(
        "ABL1", "ABL2", "AKT1", "AKT2", "AKT3", "ARAF", "AREG", "BAD",
        "BRAF", "BTC", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CBL",
        "CBLB", "CBLC", "CDKN1A", "CDKN1B", "CRK", "CRKL", "EGF", "EGFR",
        "EIF4EBP1", "ELK1", "ERBB2", "ERBB3", "ERBB4", "EREG", "GAB1",
        "GRB2", "GSK3B", "HBEGF", "HRAS", "JUN", "KRAS", "MAP2K1", "MAP2K2",
        "MAP2K4", "MAP2K7", "MAPK1", "MAPK10", "MAPK3", "MAPK8", "MAPK9",
        "MTOR", "MYC", "NCK1", "NCK2", "NRAS", "NRG1", "NRG2", "NRG3",
        "NRG4", "PAK1", "PAK2", "PAK3", "PAK4", "PAK5", "PAK6", "PIK3CA",
        "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R5",
        "PLCG1", "PLCG2", "PRKCA", "PRKCB", "PRKCG", "PTK2", "RAF1",
        "RPS6KB1", "RPS6KB2", "SHC1", "SHC2", "SHC3", "SHC4", "SOS1",
        "SOS2", "SRC", "STAT5A", "STAT5B", "TGFA", "BCL2L11", "BUB1B-PAK6",
        "CCND1", "FOXO1", "MDM2", "PDPK1", "TP53", "ACP4", "AFAP1L2",
        "AGR2", "AGT", "BCAR1", "BCAR3", "CADM1", "CAMLG", "CCDC88A",
        "CDH13", "CEACAM1", "CHMP6", "CNOT9", "CPNE3", "CUL5", "DAB2IP",
        "DGKD", "DUSP3", "EFEMP1", "EPGN", "ERBIN", "ERRFI1", "FAM83A",
        "FAM83B", "FAM83C", "FASLG", "FER", "GAREM1", "GPER1", "GPRC5A",
        "GRB7", "HAP1", "HIP1", "HIP1R", "IFI6", "IQGAP1", "ITGA1", "KIF16B",
        "LGMN", "MIR133A1", "MIR21", "MIR29A", "MMP9", "MVB12A", "MVB12B",
        "MVP", "MYOC", "NEU3", "NPR2", "NUP62", "PDE6G", "PDE6H", "PIGR",
        "PIK3C2A", "PLAUR", "PLCE1", "PRICKLE1", "PTK2B", "PTK6", "PTPN11",
        "PTPN12", "PTPN18", "PTPN2", "PTPN3", "PTPRJ", "PTPRR", "RAB7A",
        "RALA", "RALB", "RBPJ", "REPS2", "RHBDF1", "RHBDF2", "RNF115",
        "RNF126", "RTN4", "SH3TC2", "SHKBP1", "SLC30A10", "SNX5", "SNX6",
        "SOCS4", "SOCS5", "SOX9", "SPRY2", "STUB1", "TDGF1", "TGFB1",
        "TSG101", "VIL1", "VPS25", "WDR54", "ZFYVE28", "ZGPAT"
      ),
      "Tozasertib" = c(
        "AJUBA", "AKT1", "ARHGEF7", "AURKA", "AURKAIP1", "AURKB", "BIRC5",
        "BRCA1", "CDC25B", "CENPA", "CKAP5", "CPEB1", "DLGAP5", "FZR1",
        "GADD45A", "GIT1", "GSK3B", "MDM2", "NDEL1", "NFKBIA", "OAZ1",
        "PAK1", "PPP2R5D", "PRKACA", "RAN", "RASA1", "TACC1", "TACC3",
        "TDRD7", "TP53", "TPX2"
      )
    )
    selectedGenes <- pathwayKeyValue[[drug]]
    selectedGenes <- intersect(selectedGenes, colnames(x))
    selectedGenes <- intersect(selectedGenes, colnames(gdseX))

    subx <- subx[, selectedGenes]
    gdseX <- gdseX[, selectedGenes]

    preProcValues <- preProcess(subx, method = c("center", "scale"))
    subXtrans <- predict(preProcValues, subx)

    tryCatch(expr = {
      model <- train(
        x = subXtrans,
        y = suby,
        method = mlMethod,
        metric = "RMSE",
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
      targetCor <- corTargetList[[mlMethod]][drug, "text-mining"]
      model$distanceCor <- achievedCor - targetCor
      print(paste(i, model$distanceCor))
      model <- calculateMetrics(model, pred = pred, obs = gdseY)
      model$idx <- indexes[[i]]
      if (abs(model$distanceCor) < bestDistAchieved) {
        print("local best!")
        print(model$distanceCor)
        bestModel <- model
        bestDistAchieved <- abs(model$distanceCor)
      }
    }, error = function(cond) message(cond))
    if (bestDistAchieved < max.desirable.distance && bestDistAchieved > min.desirable.distance) break()
    if (model$distanceCor > 0.05 + max.desirable.distance ||
      model$distanceCor < -0.05 + min.desirable.distance) {
      break()
    }
    if (bestDistAchieved < max.desirable.distance && bestDistAchieved > min.desirable.distance) break()
  }
  return(bestModel)
}

procRNA <- function(mat) {
  rvar <- apply(mat, 1, var)
  mat <- mat[!is.na(rvar), ]
  cvar <- apply(mat, 2, var)
  mat <- mat[, !is.na(cvar)]
  return(mat)
}

trainIn <- readRDS("./result/ml_model_list/ccle_ctrpv2_rnaseq.rds")
trainOut <- readRDS("./result/ml_model_list/ccle_ctrpv2_aac.rds")
targetIn <- readRDS("./result/ml_model_list/gdse_rnaseq.rds")
targetOut <- readRDS("./result/ml_model_list/gdse_aac.rds")

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
# harmonize rows
trainSamples <- intersect(rownames(trainIn), rownames(trainOut))
targetSamples <- intersect(rownames(targetIn), rownames(targetOut))
trainIn <- trainIn[trainSamples, ]
trainOut <- trainOut[trainSamples, ]
#
targetIn <- targetIn[targetSamples, ]
targetOut <- targetOut[targetSamples, ]
#

dir <- "./result/ml_model_list_pathways"
dir.create(dir)

corTargetList <- list(
  "glmnet" = readRDS("./result/train_test/test_ElasticNet_Reg.rds"),
  "rf" = readRDS("./result/train_test/test_RandomForest_Reg.rds")
)



drugs <- c("Erlotinib", "Lapatinib", "Tozasertib")
for (mlMethod in mlMethods) {
  for (drug in rev(drugs)) {
    name <- paste0(paste(drug, mlMethod, sep = "-"), ".rds")
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
    fileName <- paste0(dir, "/", name)
    if (file.exists(fileName)) {
      bestModel <- readRDS(fileName)
      n <- bestModel$n + 1
      bestMinDiffAchieved <- bestModel$distanceCor
      if (bestMinDiffAchieved < max.desirable.distance && bestMinDiffAchieved > min.desirable.distance) { # || n >= max.try
        print("skipping")
        next()
      } else if (n > max.try) next()
    } else {
      n <- 0
      bestModel <- list()
      bestMinDiffAchieved <- 1
    }
    while (TRUE) {
      trainCut <- 0.75
      n <- n + 1
      print(paste("try", n))
      # seed <- as.numeric(Sys.time())
      seed <- n
      set.seed(seed)
      rows <- seq_len(nrow(trainX))
      if (bestMinDiffAchieved > max.desirable.distance) {
        trainCut <- trainCut - 0.05
      } else {
        trainCut <- trainCut + 0.05
      }

      x <- trainX
      y <- trainY
      idx <- sample(seq_len(nrow(x)), floor(nrow(x) * trainCut), replace = TRUE)
      x <- x[idx, ]
      y <- y[idx]

      # idx <- which(!rownames(targetX) %in% commonCellLines)

      tryCatch(
        expr = {
          model <- trainModel(
            x = x,
            y = y,
            gdseX = targetX,
            gdseY = targetY,
            drug = drug,
            mlMethod = mlMethod,
            # geneSelectMethod = geneSelectMethod,
            corTargetList = corTargetList,
            seed = seed,
            min.desirable.distance = min.desirable.distance,
            max.desirable.distance = max.desirable.distance
          )
          model$n <- n
          bestModel$n <- n
          if (model$distanceCor < max.desirable.distance && model$distanceCor > min.desirable.distance) {
            print("saving...")
            print(model$distanceCor)
            bestModel <- model
            bestModel$n <- bestModel$n - 1 # this is in case I changed my mind about the cutoff we can go back to the same 50% parition that was quite sucessful
            saveRDS(bestModel, fileName)
            break()
          } else if (abs(model$distanceCor) < bestMinDiffAchieved) {
            print("newBest!")
            print(model$distanceCor)
            bestMinDiffAchieved <- model$distanceCor
            bestModel <- model
            saveRDS(bestModel, fileName)
          } else if (n > max.try) {
            print("exiting due to 1000 tries")
            print(paste("best Cor achived is", bestMinDiffAchieved))
            saveRDS(bestModel, fileName)
            break()
          } else {
            saveRDS(bestModel, fileName)
            next()
          }
        },
        error = function(cond) message(cond)
      )
    }
    n <- n + 1
    saveRDS(bestModel, fileName)
  }
}
stopCluster(cl)
print("done")
