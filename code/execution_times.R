library(ggplot2)
library(ggpubr)
library(ggsci)
library(readxl)
df <- read_xlsx("./result/Supplementary-data-2-performance-indexes.xlsx", sheet = "execution times")
df <- df[df$"Machine learning method" == "DeepLearning", ]
df <- df[order(df$`Execution time (s) for 38 drugs`), ]
p <- ggplot(df, aes(
  x = `feature selection method`,
  y = `Execution time (s) for 38 drugs`
)) +
  geom_bar(
    stat = "identity", fill = "#E64B35FF",
    alpha = .6, width = .4
  ) +
  xlab("Feature selection methods") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic()
# theme(axis.title.x = element_text())
pdf("./result/exeuction-times-ml.pdf", height = 3, width = 4)
plot(p)
dev.off()
# "#f68060",
