library(viridis)
library(ggpubr)
library(dplyr)
library(Metrics)
library(tools)


args = commandArgs(trailingOnly=TRUE)

if (length(args) < 5) {
  stop("Please supply arguments: pSet, method, problem, metric, type, drugname (optional)", call.=FALSE)
} else if (length(args)==6) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  metric <- args[4]
  type <- args[5]
  drugname <- args[6]
} else if (length(args)==5) {
  pSet <- args[1]
  method <- args[2]
  problem <- args[3]
  metric <- args[4]
  type <- args[5]
  drugname <- NULL
}

print(sprintf("pSet: %s, method: %s, problem: %s, metric: %s, type: %s", pSet, method, problem, metric, type))

if (pSet == "CCLE"){
  trainset <- "GDSC"
  
} else if (pSet == "GDSC") {
  trainset <- "CCLE"
}

genespath <- "./genes/"

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
  tm_count <- readRDS("../data/tm_count.rds")
  
  if (problem == "class"){
    if (metric == "Accuracy"){
      tm_data <- sapply(tm, function(temp) temp$stats$overall["Accuracy"])
      top500_data <- sapply(top500, function(temp) temp$stats$overall["Accuracy"])
      top100_data <- sapply(top100, function(temp) temp$stats$overall["Accuracy"])
      ntm_data <- sapply(ntm, function(temp) temp$stats$overall["Accuracy"])
      ft_data <- sapply(ft, function(temp) temp$stats$overall["Accuracy"])
      L1000_data <- sapply(L1000, function(temp) temp$stats$overall["Accuracy"])
    }
    else if (metric == "BalancedAccuracy"){
      tm_data <- sapply(tm, function(temp) temp$stats$byClass["Balanced Accuracy"])
      top500_data <- sapply(top500, function(temp) temp$stats$byClass["Balanced Accuracy"])
      top100_data <- sapply(top100, function(temp) temp$stats$byClass["Balanced Accuracy"])
      ntm_data <- sapply(ntm, function(temp) temp$stats$byClass["Balanced Accuracy"])
      ft_data <- sapply(ft, function(temp) temp$stats$byClass["Balanced Accuracy"])
      L1000_data <- sapply(L1000, function(temp) temp$stats$byClass["Balanced Accuracy"])
    }
    else if (metric == "AUC"){
      tm_data <- sapply(tm, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
      top500_data <- sapply(top500, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
      top100_data <- sapply(top100, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
      ntm_data <- sapply(ntm, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
      ft_data <- sapply(ft, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
      L1000_data <- sapply(L1000, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
      
    }
    data <- data.frame(
      name=c(rep("500",10), rep("100",10), rep("tm",10), rep("ntm",10), rep("ft",10), rep("L1000",10)),
      value=c(unlist(top500_data), unlist(top100_data), unlist(tm_data), unlist(ntm_data), unlist(ft_data), unlist(L1000_data))
    )
    num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
    num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
    num_tm_count <- tm_count[tm_count$name == drugname,]$count
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s %s", pSet,drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s\nTM Genes: %s", num_samples_ccle, num_samples_gdsc, num_tm_count), hjust = 0, vjust = 1)
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
      name=c(rep("500",20), rep("100",20), rep("tm",20), rep("ntm",20), rep("ft",20), rep("L1000",20)),
      value=c(unlist(top500_data), unlist(top100_data), unlist(tm_data), unlist(ntm_data), unlist(ft_data), unlist(L1000_data))
    )
    num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
    num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
    num_tm_count <- tm_count[tm_count$name == drugname,]$count
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s %s", pSet,drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s\ntm genes: %s", num_samples_ccle, num_samples_gdsc, num_tm_count), hjust = 0, vjust = 1)
    return(plt)
  }
}

plotbw_test <- function(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc){
  tm <- readRDS(sprintf("../test_output/%s/%s/%s_%s_tm.rds", pSet, problem, drugname, method))
  top500 <- readRDS(sprintf("../test_output/%s/%s/%s_%s_500.rds", pSet, problem, drugname, method))
  top100 <- readRDS(sprintf("../test_output/%s/%s/%s_%s_100.rds", pSet, problem, drugname, method))
  ntm <- readRDS(sprintf("../test_output/%s/%s/%s_%s_ntm.rds", pSet, problem, drugname, method))
  ft <- readRDS(sprintf("../test_output/%s/%s/%s_%s_ft.rds", pSet, problem, drugname, method))
  L1000 <- readRDS(sprintf("../test_output/%s/%s/%s_%s_L1000.rds",pSet, problem, drugname, method))
  
  sample_count_ccle <- readRDS("../data/sample_count_ccle.rds")
  sample_count_gdsc <- readRDS("../data/sample_count_gdsc.rds")
  tm_count <- readRDS("../data/tm_count.rds")
  
  if (problem == "class"){
    if (metric == "Accuracy"){
      tm_data <- sapply(tm, function(temp) temp$stats$overall["Accuracy"])
      top500_data <- sapply(top500, function(temp) temp$stats$overall["Accuracy"])
      top100_data <- sapply(top100, function(temp) temp$stats$overall["Accuracy"])
      ntm_data <- sapply(ntm, function(temp) temp$stats$overall["Accuracy"])
      ft_data <- sapply(ft, function(temp) temp$stats$overall["Accuracy"])
      L1000_data <- sapply(L1000, function(temp) temp$stats$overall["Accuracy"])
    }
    else if (metric == "BalancedAccuracy"){
      tm_data <- sapply(tm, function(temp) temp$stats$byClass["Balanced Accuracy"])
      top500_data <- sapply(top500, function(temp) temp$stats$byClass["Balanced Accuracy"])
      top100_data <- sapply(top100, function(temp) temp$stats$byClass["Balanced Accuracy"])
      ntm_data <- sapply(ntm, function(temp) temp$stats$byClass["Balanced Accuracy"])
      ft_data <- sapply(ft, function(temp) temp$stats$byClass["Balanced Accuracy"])
      L1000_data <- sapply(L1000, function(temp) temp$stats$byClass["Balanced Accuracy"])
    }
    else if (metric == "AUC"){
      tm_data <- sapply(tm, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
      top500_data <- sapply(top500, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
      top100_data <- sapply(top100, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
      ntm_data <- sapply(ntm, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
      ft_data <- sapply(ft, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
      L1000_data <- sapply(L1000, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
      
    }
    data <- data.frame(
      name=c(rep("500",10), rep("100",10), rep("tm",10), rep("ntm",10), rep("ft",10), rep("L1000",10)),
      value=c(unlist(top500_data), unlist(top100_data), unlist(tm_data), unlist(ntm_data), unlist(ft_data), unlist(L1000_data))
    )
    num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
    num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
    num_tm_count <- tm_count[tm_count$name == drugname,]$count
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s %s", pSet,drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s\nTM Genes: %s", num_samples_ccle, num_samples_gdsc, num_tm_count), hjust = 0, vjust = 1)
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
      name=c(rep("tm",20), rep("500", 20), rep("100",20), rep("ntm",20), rep("ft", 20), rep("L1000", 20)),
      value=c(unlist(tm_data), unlist(top500_data), unlist(top100_data), unlist(ntm_data), unlist(ft_data), unlist(L1000_data))
    )
    num_samples_gdsc <- sample_count_gdsc[sample_count_gdsc$name == drugname,]$count
    num_samples_ccle <- sample_count_ccle[sample_count_ccle$name == drugname,]$count
    num_tm_count <- tm_count[tm_count$name == drugname,]$count
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6)
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("Testing on GDSC %s %s", drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", -Inf, Inf, label = sprintf("CCLE Samples: %s\nGDSC Samples: %s\nTM Genes: %s", num_samples_ccle, num_samples_gdsc, num_tm_count), hjust = 0, vjust = 1)
    return(plt)
  }
}

plot_bar <- function(pSet, drugname, metric, method, problem, sample_count_ccle, sample_count_gdsc){
  if (problem == "regression") {
    plot_data <- data.frame(Type=character(),COR=numeric()) 
    models <- list("tm", "500", "100", "ntm", "ft", "L1000")
    for (model in models){
      print(model)
      train_out <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_%s.rds", trainset, problem, drugname, method,model))
      train_data <- sapply(train_out, function(temp) temp$prediction %>% summarise(cor = list(cor(pred, obs))))
      resample <- gsub("\\..*","",names(which.max(unlist(train_data))))
      test_out <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",pSet, problem,drugname, method, model))
      test_data <- sapply(test_out, function(temp) cor(temp$pred, temp$original))
      plot_data <- plot_data %>% add_row(Type = model, COR = test_data[resample][[1]][[1]])
    }
    plt<-ggplot(data=plot_data, aes(x=Type, y=COR)) + geom_bar(stat="identity",fill="steelblue") + geom_text(aes(label=signif(COR,digits=3)), vjust=-0.3, size=3.5)
    plt<- plt + labs(title=sprintf("Testing on GDSC %s %s", drugname, metric))
    return(plt)
  } else {
    plot_data <- data.frame(Type=character(),val=numeric()) 
    models <- list("tm", "500", "100", "ntm", "ft", "L1000")
    for (model in models){
      print(model)
      if (metric == "Accuracy"){
        train_out <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_%s.rds", trainset, problem, drugname, method,model))
        train_data <- sapply(train_out, function(temp) temp$stats$overall["Accuracy"])
        resample <- names(which.max(unlist(train_data)))
        test_out <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",pSet, problem,drugname, method, model))
        test_data <- sapply(test_out, function(temp) temp$stats$overall["Accuracy"])
        plot_data <- plot_data %>% add_row(Type = model, val = test_data[resample][[1]][[1]])
      }
      else if (metric == "BalancedAccuracy"){
        train_out <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_%s.rds", trainset, problem, drugname, method,model))
        train_data <- sapply(train_out, function(temp) temp$stats$byClass["Balanced Accuracy"])
        resample <- names(which.max(unlist(train_data)))
        test_out <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",pSet, problem,drugname, method, model))
        test_data <- sapply(test_out, function(temp) temp$stats$byClass["Balanced Accuracy"])
        plot_data <- plot_data %>% add_row(Type = model, val = test_data[resample][[1]][[1]])
      }
      else if (metric == "AUC"){
        train_out <- readRDS(sprintf("../train_output/%s/%s/output/%s_%s_%s.rds", trainset, problem, drugname, method,model))
        train_data <- sapply(train_out, function(temp) temp$prediction %>% summarise(auc = list(auc(ifelse(obs == "sensitive",1,0),predict.prob.sensitive))))
        resample <- gsub("\\..*","",names(which.max(unlist(train_data))))
        test_out <- readRDS(sprintf("../test_output/%s/%s/%s_%s_%s.rds",pSet, problem,drugname, method, model))
        test_data <- sapply(test_out, function(temp) auc(ifelse(temp$pred$obs == "sensitive",1,0),temp$pred$prob$sensitive))
        plot_data <- plot_data %>% add_row(Type = model, val = test_data[resample][[1]][[1]])
      }
    }
    plt<-ggplot(data=plot_data, aes(x=Type, y=val)) + geom_bar(stat="identity",fill="steelblue") + geom_text(aes(label=signif(val,digits=3)), vjust=-0.3, size=3.5)
    plt<- plt + labs(title=sprintf("Testing on GDSC %s %s", drugname, metric))
    return(plt)
  }
}

if (is.null(drugname)) {
  files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
} else {
  files = list(drugname)
}

for (file in files){
  drugname <- file_path_sans_ext(file)
  print(drugname)
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
}
