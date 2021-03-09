library(hrbrthemes)
library(viridis)
library(ggpubr)

drugname <- "Bortezomib"
pSet <- "CCLE"
metric <- "COR"
method <- "glmnet"
problem <- "regression"

plot_univariate <- function(pSet, drugname, abs = FALSE){
  tm <- readRDS(sprintf("../train_output/univariate/%s/%s_tm.rds", pSet,drugname))
  ntm <- readRDS(sprintf("../train_output/univariate/%s/%s_ntm.rds", pSet,drugname))
  if (abs == TRUE){
    tm <- abs(tm)
    ntm <- abs(ntm)
  }
  data <- data.frame(
    name=c(rep("text_mining",length(tm)), rep("not_text_mining",length(ntm))),
    value=c(tm, ntm)
  )
  plt <- ggplot(data, aes(x=name, y=value), inherit.aes = FALSE) + geom_boxplot(alpha=0.6)
  plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s COR, abs: %s", drugname, abs),x="Feature Selection", y = "COR")
  plt <- plt + stat_compare_means(label.x = 2)
  plt <- plt + stat_compare_means(method = "t.test")
  return(plt)
}

plotbw_train <- function(pSet, drugname, metric, method, problem){
  tm <- readRDS(sprintf("../train_output/CCLE/regression/output/%s_%s_%s_%s_tm.rds", pSet, drugname, method, problem))
  top500 <- readRDS(sprintf("../train_output/CCLE/regression/output/%s_%s_%s_%s_500.rds", pSet, drugname, method, problem))
  top100 <- readRDS(sprintf("../train_output/CCLE/regression/output/%s_%s_%s_%s_100.rds", pSet, drugname, method, problem))
  ntm <- readRDS(sprintf("../train_output/CCLE/regression/output/%s_%s_%s_%s_ntm.rds", pSet,drugname, method, problem))
  ft <- readRDS(sprintf("../train_output/CCLE/regression/output/%s_%s_%s_%s_ft.rds", pSet,drugname, method, problem))
  L1000 <- readRDS(sprintf("../train_output/CCLE/regression/output/%s_%s_%s_%s_L1000.rds", pSet,drugname, method, problem))
  
  if (problem == "class"){
    if (metric == "Accuracy"){
      tm_data <- sapply(tm$stats, function(temp) temp$overall["Accuracy"])
      top500_data <- sapply(top500$stats, function(temp) temp$overall["Accuracy"])
      top100_data <- sapply(top100$stats, function(temp) temp$overall["Accuracy"])
      ntm_data <- sapply(ntm$stats, function(temp) temp$overall["Accuracy"])
    }
    else if (metric == "Balanced Accuracy"){
      tm_data <- sapply(tm$stats, function(temp) temp$byClass["Balanced Accuracy"])
      top500_data <- sapply(top500$stats, function(temp) temp$byClass["Balanced Accuracy"])
      top100_data <- sapply(top100$stats, function(temp) temp$byClass["Balanced Accuracy"])
      ntm_data <- sapply(ntm$stats, function(temp) temp$byClass["Balanced Accuracy"])
    }
    else if (metric == "AUC"){
      tm_data <- (tm$prediction %>% group_by(tm$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      top500_data <- (top500$prediction %>% group_by(top500$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      top100_data <- (top100$prediction %>% group_by(top100$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      ntm_data <- (ntm$prediction %>% group_by(ntm$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      
    }
    data <- data.frame(
      name=c(rep("text_mining",25), rep("500_genes",25), rep("100_genes",25), rep("not_text_mining",25)),
      value=c(tm_data, top500_data, top100_data, ntm_data)
    )
    
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6) 
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s", drugname, metric),x="Feature Selection", y = "Percentage")
    return(plt)
  }
  else {
    if (metric == "RMSE"){
      tm_data <- sapply(tm, function(temp) temp$stats["RMSE"])
      top500_data <- sapply(top500, function(temp) temp$stats["RMSE"])
      top100_data <- sapply(top100, function(temp) temp$stats["RMSE"])
      ntm_data <- sapply(ntm, function(temp) temp$stats["RMSE"])
      ft_data <- sapply(ft, function(temp) temp$stats["RMSE"])
      L1000_data <- sapply(L1000, function(temp) temp$stats["RMSE"])
    }
    else if (metric == "COR"){
      tm_data <- sapply(tm, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
      top500_data <- sapply(top500, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
      top100_data <- sapply(top100, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
      ntm_data <- sapply(ntm, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
      ft_data <- sapply(ft, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
      L1000_data <- sapply(L1000, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
    }
    data <- data.frame(
      name=c(rep("500_genes",25), rep("100_genes",25), rep("tm_genes",25), rep("ntm_genes",25), rep("top_cor",25), rep("L1000_genes",25)),
      value=c(unlist(top500_data), unlist(top100_data), unlist(tm_data), unlist(ntm_data), unlist(ft_data), unlist(L1000_data))
    )
    
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s %s", pSet,drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    return(plt)
  }
}

plotbw_test <- function(pSet, drugname, metric, method, problem){
  tm <- readRDS(sprintf("../test_output/GDSC2/%s_%s_%s_%s_tm.rds", pSet,drugname, method, problem))
  top500 <- readRDS(sprintf("../test_output/GDSC2/%s_%s_%s_%s_500.rds", pSet,drugname, method, problem))
  top100 <- readRDS(sprintf("../test_output/GDSC2/%s_%s_%s_%s_100.rds", pSet,drugname, method, problem))
  ntm <- readRDS(sprintf("../test_output/GDSC2/%s_%s_%s_%s_ntm.rds", pSet,drugname, method, problem))
  ft <- readRDS(sprintf("../test_output/GDSC2/%s_%s_%s_%s_ft.rds", pSet,drugname, method, problem))
  L1000 <- readRDS(sprintf("../test_output/GDSC2/%s_%s_%s_%s_L1000.rds", pSet,drugname, method, problem))
  
  if (problem == "class"){
    if (metric == "Accuracy"){
      tm_data <- sapply(tm$stats, function(temp) temp$overall["Accuracy"])
      top500_data <- sapply(top500$stats, function(temp) temp$overall["Accuracy"])
      top100_data <- sapply(top100$stats, function(temp) temp$overall["Accuracy"])
      ntm_data <- sapply(ntm$stats, function(temp) temp$overall["Accuracy"])
    }
    else if (metric == "Balanced Accuracy"){
      tm_data <- sapply(tm$stats, function(temp) temp$byClass["Balanced Accuracy"])
      top500_data <- sapply(top500$stats, function(temp) temp$byClass["Balanced Accuracy"])
      top100_data <- sapply(top100$stats, function(temp) temp$byClass["Balanced Accuracy"])
      ntm_data <- sapply(ntm$stats, function(temp) temp$byClass["Balanced Accuracy"])
    }
    else if (metric == "AUC"){
      tm_data <- (tm$prediction %>% group_by(tm$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      top500_data <- (top500$prediction %>% group_by(top500$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      top100_data <- (top100$prediction %>% group_by(top100$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      ntm_data <- (ntm$prediction %>% group_by(ntm$prediction$resample) %>% summarise(auc = as.double(auc(original,predict.prob.resistant))))$auc
      
    }
    data <- data.frame(
      name=c(rep("text_mining",25), rep("500_genes",25), rep("100_genes",25), rep("not_text_mining",25)),
      value=c(tm_data, top500_data, top100_data, ntm_data)
    )
    
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6) 
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s", drugname, metric),x="Feature Selection", y = "Percentage")
    return(plt)
  }
  else {
    if (metric == "RMSE"){
      tm_data <- sapply(tm, function(temp) temp$stats["RMSE"])
      top500_data <- sapply(top500, function(temp) temp$stats["RMSE"])
      top100_data <- sapply(top100, function(temp) temp$stats["RMSE"])
      ntm_data <- sapply(ntm, function(temp) temp$stats["RMSE"])
      ft_data <- sapply(ft, function(temp) temp$stats["RMSE"])
      L1000_data <- sapply(L1000, function(temp) temp$stats["RMSE"])
    }
    else if (metric == "COR"){
      tm_data <- (tm$prediction %>% group_by(tm$prediction$resample) %>% summarise(cor = cor(pred, obs)))$cor
      top500_data <- (top500$prediction %>% group_by(top500$prediction$resample) %>% summarise(cor = cor(pred, obs)))$cor
      top100_data <- (top100$prediction %>% group_by(top100$prediction$resample) %>% summarise(cor = cor(pred, obs)))$cor
      ntm_data <- (ntm$prediction %>% group_by(ntm$prediction$resample) %>% summarise(cor = cor(pred, obs)))$cor
      ft_data <- (ft$prediction %>% group_by(ft$prediction$resample) %>% summarise(cor = cor(pred, obs)))$cor
      L1000_data <- (L1000$prediction %>% group_by(L1000$prediction$resample) %>% summarise(cor = cor(pred, obs)))$cor
    }
    data <- data.frame(
      name=c(rep("text_mining",25), rep("500_genes", 25), rep("100_genes",25), rep("not_text_mining",25), rep("top_cor", 25), rep("L1000", 25)),
      value=c(tm_data, top500_data, top100_data, ntm_data, ft_data, L1000_data)
    )
    
    print(tm_data)
    
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s", drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    return(plt)
  }
}

plotbw_100 <- function(pSet, drugname, metric, method, problem){
  tm_ft <- readRDS(sprintf("../train_output/%s_%s_%s_%s_tm_ft.rds", pSet,drugname, method, problem))
  ft_100 <- readRDS(sprintf("../train_output/%s_%s_%s_%s_ft_100.rds", pSet,drugname, method, problem))
  tm_data <- sapply(tm_ft$output, function(temp) temp$stats["RMSE"])
  ft_data <- sapply(ft_100$output, function(temp) temp$stats["RMSE"])
  data <- data.frame(
    name=c(rep("text_mining",10), rep("top_cor", 10)),
    value=c(tm_data, ft_data)
  )
  plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
  plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s", drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
  return(plt)
}

drugname <- "Alisertib"
plot<- plotbw_train(pSet, drugname, metric, method, problem)
print(plot)

pdf(sprintf("../result/test.pdf", pSet,drugname, metric, problem))
print(plot)
dev.off()
