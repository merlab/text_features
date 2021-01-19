drugname <- "Nilotinib"
pSet <- "GDSC2"
metric <- "AUC"
method <- "glmnet"
problem <- "class"

plotbw <- function(pSet, drugname, metric, method, problem){
  tm <- readRDS(sprintf("../train_output/%s_%s_%s_%s_tm.rds", pSet,drugname, method, problem))
  top500 <- readRDS(sprintf("../train_output/%s_%s_%s_%s_500.rds", pSet,drugname, method, problem))
  top100 <- readRDS(sprintf("../train_output/%s_%s_%s_%s_100.rds", pSet,drugname, method, problem))
  ntm <- readRDS(sprintf("../train_output/%s_%s_%s_%s_ntm.rds", pSet,drugname, method, problem))
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
    tm_data <- sapply(tm$stats, function(temp) temp["RMSE"])
    top500_data <- sapply(top500$stats, function(temp) temp["RMSE"])
    top100_data <- sapply(top100$stats, function(temp) temp["RMSE"])
    ntm_data <- sapply(ntm$stats, function(temp) temp["RMSE"])
    
    tmcor <- cor(tm$prediction$pred, tm$prediction$obs)
    top500cor <- cor(top500$prediction$pred, top500$prediction$obs)
    top100cor <- cor(top100$prediction$pred, top100$prediction$obs)
    ntmcor <- cor(ntm$prediction$pred, ntm$prediction$obs)
    
    cordata <- data.frame(
      name=c("text_mining", "500_genes", "100_genes", "not_text_mining"),
      value=c(round(tmcor, digits = 3), round(top500cor, digits = 3), round(top100cor, digits=3), round(ntmcor,digits=3)))
    
    data <- data.frame(
      name=c(rep("text_mining",25), rep("500_genes",25), rep("100_genes",25), rep("not_text_mining",25)),
      value=c(tm_data, top500_data, top100_data, ntm_data)
    )
    
    plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6) +  geom_text(data = cordata, aes(label=value, y = max(data["value"]) + 0.01))
    plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s", drugname, metric),x="Feature Selection", y = sprintf("%s", metric))
    plt <- plt + annotate("text", label = "COR values:",y=max(data["value"] + 0.016), x = "100_genes")
    return(plt)
  }
}

plot <- plotbw(pSet,drugname,metric, method, problem)
print(plot)
pdf(sprintf("../result/%s_%s_%s_%s.pdf", pSet,drugname, metric, problem))
print(plot)
dev.off()
