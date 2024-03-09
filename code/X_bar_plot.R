library(ggplot2)
library(readxl)
library(ggpubr)
library(ggsci)
fPath <- "./mlModelMetrics_pathway.xlsx"
fTxt <- "./result/Supplementary-data-2-performance-indexes.xlsx"
sheets <- c("ElasticNet perf metrics-test", "RandomForest perf metrics-test")
p <- list()
for (i in sheets) {
  path <- read_xlsx(fPath, sheet = i)[, c("drug", "pearsonCor", "feature selection method")]
  txt <- read_xlsx(fTxt, sheet = i)
  txt <- txt[
    txt$drug %in% path$drug & txt$"feature selection method" == "text-mining",
    c("drug", "pearsonCor", "feature selection method")
  ]
  d <- rbind(path, txt)
  d$`feature selection method` <- factor(d$`feature selection method`, levels = c("text-mining", "pathway"))
  p[[i]] <- ggplot(d, aes(x = drug, y = pearsonCor, fill = `feature selection method`)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.75) +
    scale_fill_npg(
      labels = c("Text-mining", "Pathway"),
      name = "Feature selection method"
    ) +
    ggtitle(gsub("perf metrics-test", "", i)) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    xlab("Drug") +
    ylab("Pearson correlation") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic()
}
pdf("./ML_results_path.pdf", height = 3, width = 6)
plot(ggarrange(plotlist = p, common.legend = TRUE, nrow = 1, ncol = 2))
dev.off()
