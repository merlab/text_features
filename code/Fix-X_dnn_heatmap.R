library(readxl)
library(ggplot2)
library(ggpubr)
library(writexl)
library(tidyr)
library(ComplexHeatmap)
library(viridis)

formatXLSX <- function(x, valCol = "pearsonCor") {
  x <- x[, c("drug", "feature selection method", valCol)] # "feature selection method"
  colnames(x) <- c("drug", "Type", "value")
  y <- tidyr::pivot_wider(x, names_from = "Type", values_from = "value")
  y <- as.data.frame(y)
  rownames(y) <- y$drug
  y$drug <- NULL
  return(y)
}
h <- 5
methods <- c("spearmanCor", "kendallCor", "RMSE", "MSE", "MAE") #
f <- "./result/Supplementary-data-2-performance-indexes.xlsx"
pdf("./DNN_results_heat.pdf", width = h * 1.4 * 2, height = h)
for (i in methods) {
  dnn <- read_xlsx(f, sheet = "DeepLearning perf metrics")
  dnn <- formatXLSX(dnn, valCol = i)
  dnn <- data.matrix(dnn)
  dnn <- (dnn - mean(dnn)) / sd(dnn)
  plot(Heatmap(t(dnn),
    cluster_columns = TRUE,
    cluster_rows = FALSE,
    row_title = i
  ))
}
dev.off()
df <- data.frame()
for (i in methods) {
  dnn <- read_xlsx(f, sheet = "DeepLearning perf metrics")
  dnn <- formatXLSX(dnn, valCol = i)
  dnn <- data.matrix(dnn)
  for (j in colnames(dnn)) {
    if (j == "text-mining") next()
    p <- wilcox.test(x = dnn[, "text-mining"], y = dnn[, j])$p.value
    df <- rbind(df, c(i, j, p))
  }
}
colnames(df) <- c("feature", "metric", "pval")
df[, "pval"] <- as.numeric(df[, "pval"])
df$log10pval <- -log10(df[, "pval"])
# df <- df[df[, "pval"] <= 0.05, ]

p <- ggplot(df, aes(x = feature, y = metric, size = log10pval, color = log10pval)) +
  geom_point() +
  scale_color_viridis(
    name = "P value",
    # limits = c(.1, .4),
    # breaks = seq(.1, .4, .1)
  ) +
  scale_size(
    name = "P value",
    breaks = c(1, 1.3, 2, 4, 8),
    labels = c(
      "0.1",
      "0.05",
      expression(1 * "\u00D7" * 10^"-2"),
      expression(1 * "\u00D7" * 10^"-4"),
      expression("<" * 1 * "\u00D7" * 10^"-8")
    ),
    range = c(3, 8.5)
  ) +
  theme_bw() +
  # coord_flip() +
  # scale_x_continuous(
  #   position = "bottom",
  #   expand = c(0, 0),
  #   breaks = seq(1, max(df$x)),
  #   limits = c(0.5, max(df$x) + 0.5),
  #   labels = levels(df$metricClean),
  #   minor_breaks = seq(0.5, max(df$x))
  # ) +
  # scale_y_continuous(
  #   position = "right",
  #   expand = c(0, 0),
  #   breaks = seq(1, max(df$y)),
  #   limits = c(0.5, max(df$y) + 0.5),
  #   labels = levels(df$drugCancerCombo),
  #   minor_breaks = seq(0.5, max(df$y))
  # ) +
  xlab("") +
  ylab("")
p <- p + theme(
  axis.text.x = element_text(
    angle = 45, # 45
    hjust = -0.05,
    color = "black",
    size = 10
  ),
  axis.text.y = element_text(color = "black", size = 10),
  panel.grid.major = element_blank(),
  axis.ticks = element_blank(),
  panel.grid.minor = element_line(color = "white", linewidth = 1),
  panel.background = element_rect(fill = "#E2E2E2", color = "#E2E2E2")
)
pdf("./DNN_results_heat2.pdf", width = h * 1.4, height = h)
plot(p)
dev.off()
print("done")
