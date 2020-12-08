plotbw <- function(name, bwdata1, bwdata2, bwdata3, bwdata4){
  data <- data.frame(
    name=c(rep("text_mining",10), rep("500_genes",10), rep("100_genes",10), rep("not_text_mining",10)),
    value=c(bwdata1, bwdata2, bwdata3, bwdata4)
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
