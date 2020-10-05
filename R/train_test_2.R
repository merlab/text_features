library(caret)
library(PharmacoGx)
source("summarizeData.R")
GDSC2 <- readRDS("../data/GDSC2.rds")
df = summarizeData(pSet=GDSC2, mDataType = "Kallisto_0.46.1.rnaseq", drug = "Erlotinib", sensitivity.measure="aac_recomputed")
df = df[1:5000,]
NAsamples <- apply(assay(df), 2, function(i)any(is.na(i)))
df <- df[, !NAsamples]
df = df[, !is.na(colData(df)$aac_recomputed)]
gene_vars <- apply(assay(df), 1, var)
df <- df[gene_vars>0, ]

x <- t(assay(df))
y <- colData(df)$aac_recomputed


# exprs <- assays(df)@listData$exprs
# labels <- colData(df)$aac_recomputed
# labels <- as.data.frame(t(as.matrix(labels)))
# # removing rows with NA values from expression data 
# naSample <- apply(exprs, 2, function(x) all(is.na(x)))
# exprs <- exprs[, !naSample]
# # removing NA values from auc
# naLabel <- apply(labels, 2, function(x) all(is.na(x)))
# exprs <- exprs[,!naLabel]
# labels <- labels[,!naLabel]
# #removing features with zero variance
# gene_vars <- apply(exprs, 2, var)
# gene_names <- names(gene_vars[gene_vars > 0])
# exprs <- exprs[, gene_names]
# exprs <- as.data.frame(t(as.matrix(exprs)))
# labels <- as.data.frame(t(as.matrix(labels)))
train_result_1 <- train(x=exprs,
                        y=labels,
                        method="glmnet",
                        preProcess=c("center", "scale"),
                        maximize = TRUE,
                        na.rm = TRUE,
                        tuneGrid=expand.grid(alpha=seq(0, 1, 0.2),
                                             lambda=seq(0, 10, 1)),
                        trControl=trainControl(method="repeatedcv",
                                               number=10,
                                               repeats=5,
                                               search="grid",
                                               predictionBounds = c(0, 1),
                                               allowParallel = TRUE,
                                               verboseIter=TRUE))
GDSC2_pred_df <- data.frame(predict=predict(train_result, GDSC2_data), original=GDSC2_labels)
