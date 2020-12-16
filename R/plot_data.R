plotbw <- function(name){
  tm <- readRDS("../train_output/model_Lapatinib_tm.rds")
  top500 <- readRDS("../train_output/model_Lapatinib_500.rds")
  top100 <- readRDS("../train_output/model_Lapatinib_100.rds")
  ntm <- readRDS("../train_output/model_Lapatinib_ntm.rds")
  tm_data <- sapply(tm$stats, function(temp) temp$overall["Accuracy"])
  top500_data <- sapply(top500$stats, function(temp) temp$overall["Accuracy"])
  top100_data <- sapply(top100$stats, function(temp) temp$overall["Accuracy"])
  ntm_data <- sapply(ntm$stats, function(temp) temp$overall["Accuracy"])
  data <- data.frame(
    name=c(rep("text_mining",10), rep("500_genes",10), rep("100_genes",10), rep("not_text_mining",10)),
    value=c(tm_data, top500_data, top100_data, ntm_data)
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
    ggtitle(sprintf("%s - Training Results", name)) +
    xlab("")
}
