drugname <- "Dasatinib"
pSet <- "GDSC2"
metric <- "Balanced Accuracy"
method <- "glmnet"

plotbw <- function(pSet, drugname, metric, method){
  tm <- readRDS(sprintf("../train_output/%s_%s_%s_tm.rds", pSet,drugname, method))
  print(tm$metadata$label)
  top500 <- readRDS(sprintf("../train_output/%s_%s_%s_500.rds", pSet,drugname, method))
  top100 <- readRDS(sprintf("../train_output/%s_%s_%s_100.rds", pSet,drugname, method))
  ntm <- readRDS(sprintf("../train_output/%s_%s_%s_ntm.rds", pSet,drugname, method))
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
  data <- data.frame(
    name=c(rep("text_mining",10), rep("500_genes",10), rep("100_genes",10), rep("not_text_mining",10)),
    value=c(tm_data, top500_data, top100_data, ntm_data)
  )
  
  plt <- ggplot(data, aes(x=name, y=value, fill=name)) + geom_boxplot(alpha=0.6) 
  plt <- plt + theme(legend.position="none") + labs(title=sprintf("%s %s", drugname, metric),x="Feature Selection", y = "Percentage")
  return(plt)
}

plot <- plotbw(pSet,drugname,metric, method)
print(plot)
pdf(sprintf("../result/%s_%s_%s.pdf", pSet,drugname, metric))
print(plot)
dev.off()
