library(caret);library(PharmacoGx)
GDSC2 <- readRDS("../data/GDSC2.rds")
molprof <- summarizeMolecularProfiles(GDSC2,
                                      mDataType="Kallisto_0.46.1.rnaseq.counts",
                                      fill.missing=TRUE,
                                      summarize=TRUE,
                                      verbose=FALSE)
naSample <- apply(assay(molprof), 2, function(x) all(is.na(x)))
molprof <- molprof[, !naSample]
auc <- summarizeSensitivityProfiles(GDSC2,
                                    sensitivity.measure="aac_recomputed",
                                    fill.missing=TRUE,
                                    verbose=FALSE)
sampletotake <- intersect(rownames(colData(molprof)), dimnames(auc)[[2]])
molprof <- molprof[, sampletotake]
auc <- auc[, sampletotake]
GDSC2 <- list("molprof"=molprof, "auc"=auc)
GDSC2_drug_ind <- match("Erlotinib", dimnames(GDSC2[["auc"]])[[1]])
labels <- GDSC2$auc[GDSC2_drug_ind, ][!is.na(GDSC2$auc[GDSC2_drug_ind, ])]
features <- assay(molprof[, names(labels)])
features <- as.data.frame(t(features))
features$labels <- labels[rownames(features)]
GDSC2_data <- features
GDSC2_labels <- GDSC2_data$labels
names(GDSC2_labels) <- dimnames(GDSC2_data)[[1]]
GDSC2_data$labels <- NULL
gene_vars <- apply(GDSC2_data, 2, var)
gene_names <- names(gene_vars[gene_vars > 0])
GDSC2_data <- GDSC2_data[, gene_names]
GDSC2_data <- as.matrix(GDSC2_data[, gene_names])
x <- GDSC2_data[,1:100]
y <- GDSC2_labels
train_result_sample <- train(x=x,
                      y=y,
                      method="glmnet",
                      preProcess=c("center", "scale"),
                      metric="CI",
                      maximize=TRUE,
                      tuneGrid=expand.grid(alpha=seq(0, 1, 0.2),
                                           lambda=seq(0, 10, 1)),
                      trControl=trainControl(method="repeatedcv",
                                             number=10,
                                             repeats=5,
                                             search="grid",
                                             predictionBounds = c(0, 1),
                                             allowParallel = FALSE,
                                             verboseIter=TRUE))
GDSC2_pred_df_allfeatures <- data.frame(predict=predict(train_result_allfeatures, GDSC2_data), original=GDSC2_labels)
