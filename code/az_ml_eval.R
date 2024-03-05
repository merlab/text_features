library(ggsci)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

dir <- "./result/ml_model_list"
fs <- list.files(dir, pattern = ".*-.*.rds")
distV <- c()
for (i in fs) {
  distV[i] <- readRDS(sprintf("%s/%s", dir, i))$distanceCor
}
ns <- fs
ns <- gsub("\\.rds", "", ns)


mlModelV <- ifelse(grepl("glmnet", ns), "glmnet", "rf")
drugV <- sapply(strsplit(ns, "-"), function(x) {
  return(x[1])
})
selectMethodV <- sapply(strsplit(ns, "-"), function(x) {
  return(paste(x[-c(1, 2)], collapse = "-"))
})
plotdf <- data.frame(
  dist = distV,
  mlModel = mlModelV,
  drug = drugV,
  selectMethod = selectMethodV
)

p <- ggplot(data = plotdf, aes(x = dist, group = selectMethod, fill = selectMethod)) +
  geom_density(adjust = 1.5, alpha = .75) +
  theme_classic() +
  scale_color_npg() +
  scale_fill_npg() +
  geom_vline(xintercept = 0.01, color = "red", linetype = "dashed") +
  geom_vline(xintercept = -0.01, color = "red", linetype = "dashed") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~ mlModel + selectMethod)
#
pdf("res.pdf", height = 10, width = 10)
plot(p)
dev.off()
print("done")
