library(viridis)
library(ggpubr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 6) {
  stop("Please supply arguments: pSet, method, problem, drugname, metric, type", call.=FALSE)
} else if (length(args)==6) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  drugname <- args[4]
  metric <- args[5]
  type <- args[6]
}

print(sprintf("pSet: %s, method: %s, problem: %s, drugname: %s, metric: %s, type: %s", pSet, method, problem, drugname, metric, type))

if (pSet == "CCLE"){
  trainset <- "GDSC"
  
} else if (pSet == "GDSC") {
  trainset <- "CCLE"
}

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

plotbw_train <- function(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc){
  tm <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_tm.rds", pSet, problem, drugname, method))
  top500 <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_500.rds", pSet, problem, drugname, method))
  top100 <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_100.rds", pSet, problem, drugname, method))
  ntm <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_ntm.rds", pSet, problem, drugname, method))
  ft <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_ft.rds", pSet, problem, drugname, method))
  L1000 <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_L1000.rds",pSet, problem, drugname, method))
  
  sample_count_ccle <- readRDS("../data/sample_count_ccle.rds")
  sample_count_gdsc <- readRDS("../data/sample_count_gdsc.rds")
  
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
    num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
    num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s %s", pSet,drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s", num_samples_ccle, num_samples_gdsc), hjust = 0, vjust = 1)
    return(plt)
  }
}

plotbw_test <- function(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc){
  tm <- readRDS(sprintf("../test_output/%s/%s/%s_%s_tm.rds", pSet, problem, drugname, method))
  top500 <- readRDS(sprintf("../test_output/%s/%s/%s_%s_500.rds",pSet, problem, drugname, method))
  top100 <- readRDS(sprintf("../test_output/%s/%s/%s_%s_100.rds", pSet, problem, drugname, method))
  ntm <- readRDS(sprintf("../test_output/%s/%s/%s_%s_ntm.rds",pSet, problem, drugname, method))
  ft <- readRDS(sprintf("../test_output/%s/%s/%s_%s_ft.rds",pSet, problem, drugname, method))
  L1000 <- readRDS(sprintf("../test_output/%s/%s/%s_%s_L1000.rds", pSet, problem, drugname, method))
  
  sample_count_ccle <- readRDS("../data/sample_count_ccle.rds")
  sample_count_gdsc <- readRDS("../data/sample_count_gdsc.rds")
  
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
      tm_data <- sapply(tm, function(temp) cor(temp$pred, temp$original))
      top500_data <- sapply(top500, function(temp) cor(temp$pred, temp$original))
      top100_data <- sapply(top100, function(temp) cor(temp$pred, temp$original))
      ntm_data <- sapply(ntm, function(temp) cor(temp$pred, temp$original))
      ft_data <- sapply(ft, function(temp) cor(temp$pred, temp$original))
      L1000_data <- sapply(L1000, function(temp) cor(temp$pred, temp$original))
    }
    data <- data.frame(
      name=c(rep("text_mining",25), rep("500_genes", 25), rep("100_genes",25), rep("not_text_mining",25), rep("top_cor", 25), rep("L1000", 25)),
      value=c(unlist(tm_data), unlist(top500_data), unlist(top100_data), unlist(ntm_data), unlist(ft_data), unlist(L1000_data))
    )
    num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
    num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("Testing on GDSC %s %s", drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s", num_samples_ccle, num_samples_gdsc), hjust = 0, vjust = 1)
    return(plt)
  }
}

plot_bar <- function(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc){
  plot_data <- data.frame(Type=character(),COR=numeric()) 
  models <- list("tm", "500", "100", "ntm", "ft", "L1000")
  for (model in models){
    print(model)
    train_out <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_%s.rds", trainset, problem, drugname, method,model))
    train_data <- sapply(train_out, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
    resample <- names(which.max(unlist(train_data)))
    test_out <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",pSet, problem,drugname, method, model))
    train_data <- sapply(train_out, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
    plot_data <- plot_data %>% add_row(Type = model, COR = train_data[resample][[1]][[1]])
  }
  plt<-ggplot(data=plot_data, aes(x=Type, y=COR)) + geom_bar(stat="identity",fill="steelblue") + geom_text(aes(label=signif(COR,digits=3)), vjust=-0.3, size=3.5)
  plt<- plt + labs(title=sprintf("Testing on GDSC %s %s", drugname, metric))
  return(plt)
}

plotbw_var <- function(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc){
  top50 <- readRDS(sprintf("../train_output/misc/output/%s_%s_%s_50.rds", drugname, method, problem))
  top100 <- readRDS(sprintf("../train_output/misc/model/%s_%s_%s_100_fix.rds", drugname, method, problem))
  top200 <- readRDS(sprintf("../train_output/misc/model/%s_%s_%s_200_fix.rds", drugname, method, problem))
  top500 <- readRDS(sprintf("../train_output/misc/output/%s_%s_%s_500_fix.rds", drugname, method, problem))
  top1000 <- readRDS(sprintf("../train_output/misc/output/%s_%s_%s_1000_fix.rds", drugname, method, problem))
  
  top50_data <- sapply(top50, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
  top100_data <- sapply(top100, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
  top200_data <- sapply(top200, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
  top500_data <- sapply(top500, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
  top1000_data <- sapply(top1000, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))

  data <- data.frame(
    name=c(rep("50_genes",25), rep("100_genes", 25), rep("200_genes",25), rep("500_genes",25), rep("1000_genes", 25)),
    value=c(unlist(top50_data), unlist(top100_data), unlist(top200_data), unlist(top500_data), unlist(top1000_data))
  )
  level_order <- c("50_genes", "100_genes","200_genes","500_genes","1000_genes" )
  num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
  num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
  plt <- ggplot(data, aes(x=factor(name, levels = level_order), y=value, fill=name)) + geom_boxplot(alpha=0.6)
  plt <- plt + theme(legend.position="none") + labs(title=sprintf("Training with %s %s", drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
  plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s", num_samples_ccle, num_samples_gdsc), hjust = 0, vjust = 1)
  return(plt)
  
}

if (type == "bar"){
  plot<- plot_bar(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc)
}else if (type == "train"){
  plot<- plotbw_train(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc)
} else if (type == "test"){
  plot<- plotbw_test(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc)
} 
pdf(sprintf("../result/%s/%s/%s/%s_%s_%s.pdf", type, pSet, problem, drugname, method, metric))
print(plot)
dev.off()
