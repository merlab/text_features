library(caret)
library(psych)
library(PharmacoGx)
library(SummarizedExperiment)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(sqldf)
library(data.table)
source("./summarizeData.R")
GDSC2 <- readRDS("../data/GDSC2.rds")
drugname <- "erlotinib"

generate_df <- function(pSet, mDataType, drug){
  #create df
  df = summarizeData(pSet=pSet, mDataType = mDataType, drug = drug, sensitivity.measure="aac_recomputed")
  #df = summarizeData(pSet=GDSC2, mDataType = "Kallisto_0.46.1.rnaseq", drug = "Erlotinib", sensitivity.measure="aac_recomputed")
  df = df[, !is.na(colData(df)$aac_recomputed)]
  
  #remove samples with NA value
  NAsamples <- apply(assay(df), 2, function(i) any(is.na(i)))
  df <- df[, !NAsamples]
  
  #remove samples with no variance
  gene_vars <- apply(assay(df), 1, var)
  gene_vars <- sort(gene_vars, decreasing = TRUE)
  gene_vars <- gene_vars[gene_vars>0]
  df <- df[names(gene_vars), ]
  
  return(df)
}

subset_by_feat <- function(df, drug, textmining, subset_size){
  
  #select mined genes
  if (textmining == TRUE){
    minedgenes <- readRDS(sprintf("./%s.rds",toupper(drug)))
    dfout = df[rowData(df)$gene_name %in% minedgenes$Symbol, ]
  } else {
    dfout= df[!(rowData(df)$gene_name %in% minedgenes$Symbol), ]
    dfout = dfout[1:subset_size]
  }
  
  #produce x and y, turn y to discrete
  x <- t(assay(dfout))
  y <- colData(dfout)$aac_recomputed
  y <- ifelse(y >= 0.1,"sensitive","resistant")
  
  #print class counts for y
  counts <- table(y)
  sprintf("Resistant: %d", counts[1])
  sprintf("Sensitive: %d", counts[2])
  
  output <- list("X"=x, "Y"=y)
  return(output)
}

#second function trains model based on x and y
trainmodel <- function(x,y,name){
  trainIndex <- createDataPartition(y, p = .8, 
                                    list = TRUE, 
                                    times = 10)
  
  tgrid <- expand.grid(alpha=seq(0, 1, 0.2),
                       lambda=seq(0, 10, 1))
  
  tcontrol <- trainControl(method="repeatedcv",
                           number= 4,
                           repeats = 4,
                           search="grid",
                           savePredictions ="all",
                           allowParallel = TRUE,
                           verboseIter=TRUE,
                           classProbs = TRUE)
  
  pred_sample <- data.frame()
  result <- data.frame()
  for(res in names(trainIndex))
  {
    trIndx <- trainIndex[[res]]
    tsIndx <- setdiff(1:nrow(x), trIndx)
    
    train_result_sample <- train(x=x[trIndx, ], y=y[trIndx],
                                 method="glmnet",
                                 #preProcess=c("center", "scale"),
                                 maximize = TRUE,
                                 tuneGrid=tgrid,
                                 trControl=tcontrol)
    train_result_sample$trainingData <- NULL
    pred_sample_n <- data.frame(index = tsIndx,
                                predict=predict(train_result_sample, x[tsIndx, ]),
                                original=y[tsIndx], 
                                resample=res)
    acc <- sum(pred_sample_n$original == pred_sample_n$predict)/nrow(pred_sample_n)
    result <- rbind(result, data.frame(ACC=acc, resample=res))
    pred_sample <- rbind(pred_sample, pred_sample_n)
  }
  
  #saveRDS(train_result_sample, sprintf("model_%s.rds", name))
  #saveRDS(pred_sample, sprintf("pred_result_%s.rds", name))
  metadata <- list(table(y), dim(x))
<<<<<<< HEAD
  output <- list("pred_sample"=pred_sample, metadata)
=======
  output <- list(train_result_sample, pred_sample, metadata, result)
>>>>>>> e0bb476d209a891252668071192b71de328ef09d
  return(output)
}

df <- generate_df(GDSC2, "Kallisto_0.46.1.rnaseq", drugname)

# for text mining genes
temp1 <- subset_by_feat(df, drugname, TRUE, 500)
result1 <- trainmodel(temp1$X, temp1$Y, "erlotinib")
DT1 <- data.table(result1$pred_sample)
bwdata1 <- sqldf('SELECT count(*)*100/(SELECT count(*) FROM DT1 GROUP BY resample having predict = original) as accuracy FROM DT1 where predict = original GROUP BY resample')

# for top 500 genes
temp2 <- subset_by_feat(df, "Erlotinib", FALSE, 500)
result2 <- trainmodel(temp2$X, temp2$Y, "erlotinib")
DT2 <- data.table(result2$pred_sample)
bwdata2 <- sqldf('SELECT count(*)*100/(SELECT count(*) FROM DT2 GROUP BY resample having predict = original) as accuracy FROM DT2 where predict = original GROUP BY resample')

# for top 100 genes
temp3 <- subset_by_feat(df, "Erlotinib", FALSE, 100)
result3 <- trainmodel(temp3$X, temp3$Y, "erlotinib")
DT3 <- data.table(result3$pred_sample)
bwdata3 <- sqldf('SELECT count(*)*100/(SELECT count(*) FROM DT3 GROUP BY resample having predict = original) as accuracy FROM DT3 where predict = original GROUP BY resample')

# generate bwplot
data <- data.frame(
  name=c( rep("text_mining",10), rep("500_genes",10), rep("100_genes",10)),
  value=c( bwdata1$accuracy, bwdata2$accuracy, bwdata3$accuracy )
)

data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme(text=element_text(size=16,  family="serif")) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Training Results") +
  xlab("")



####
minedgenes <- readRDS("./ERLOTINIB.rds")
dftext = df[rowData(df)$gene_name %in% minedgenes$Symbol, ]
dfOther= df[!(rowData(df)$gene_name %in% minedgenes$Symbol), ]


