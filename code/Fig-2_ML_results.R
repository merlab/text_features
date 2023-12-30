# NOTE: this figure was modified from the old script to include the new
# results (REF, MRMR, and GA)
library(readxl)
library(ggplot2)
library(ggpubr)
# library(gridExtra)
# library(readr); library(tidyr)
# library(tools); library(Hmisc) ; library(plyr); library(RColorBrewer)
# library(reshape2); ; library(ggrepel)

# source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
source("./code/geom_flat_violin.R")

raincloud_theme <- theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5, color = "black"),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 16),
  legend.position = "right",
  plot.title = element_text(lineheight = .8, face = "bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = "black", size = 0.5, linetype = "solid"),
  axis.line.y = element_line(colour = "black", size = 0.5, linetype = "solid")
)

## -------------------------------------------------------------

create_plot <- function(dft, title = NA) {
  # fetMeth <- colnames(df)[colnames(df) != "drug"]
  dft <- as.data.frame(dft)
  #   fetMeth <- c("var-100", "var-500", "L1000-tm", "L1000", "cor-500", "text-mining")
  fetMeth <- c(
    "var-500", "var-100", "L1000", "L1000-tm", "cor-500",
    "RFE", "MRMR", "GA",
    "text-mining"
  )

  dft$Type <- factor(dft$Type, levels = fetMeth)

  # cl <- c("#B09C85FF", "#91D1C2FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF", "#DC0000FF")
  # cl <- c(
  #   "#B09C85FF", "#91D1C2FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF",
  #   "#DC0000FF", "#173dd3", "#29df49", "#24dd9d"
  # )
  # cl <- c(
  #   "#B09C85FF", "#91D1C2FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF",
  #   "#DC0000FF", "#85469d", "#bf4aa5", "#446e35"
  # )
  cl <- c(
    "#B09C85FF", "#91D1C2FF", "#4DBBD5FF", "#3C5488FF", "#00A087FF",
    "#85469d", "#bf4aa5", "#446e35", "#DC0000FF"
  )
  names(cl) <- fetMeth

  my_comparisons <- lapply(
    rev(fetMeth[1:(length(fetMeth) - 1)]),
    function(i) c(i, fetMeth[length(fetMeth)])
  )

  # df$drug <- rownames(df)

  # dft <- reshape2::melt(df, id.vars = "drug", variable.name = "Type")
  # print(head(dft))
  # stop()

  g <- ggplot(data = dft, aes(y = value, x = Type, fill = Type)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .95) +
    geom_point(aes(y = value, color = Type),
      position = position_jitter(width = .15), size = .5, alpha = 0.45
    ) +
    geom_boxplot(
      width = .1, # guides = "none",
      outlier.shape = NA, alpha = 0.5
    ) +
    expand_limits(x = 5.25) +
    guides(fill = "none") +
    guides(color = "none") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    raincloud_theme +
    xlab("") + # xlab("Feature Selection Method") +
    ylab("Pearson correlation")

  g <- g + scale_color_manual(values = cl) + scale_fill_manual(values = cl)
  g <- g + stat_compare_means(
    comparisons = my_comparisons, method = "wilcox.test", paired = TRUE,
    # label.y = 0.7,
    size = 2.5
  )
  g <- g + scale_y_continuous(breaks = seq(-0.4, 1, 0.2))
  g + theme(
    axis.text.x = element_text(
      face = "bold", # color="#993333",
      color = "black",
      size = 8, angle = -0
    ),
    axis.text.y = element_text(
      face = "bold", # color="#993333",
      color = "black",
      size = 9, angle = 0
    ),
    axis.title.y = element_text(face = "bold", size = 10)
  )
  return(g)
}

## --------------------------------------------------------------
## ---- Random-forest training and test plot --------------------

trRF <- read_xlsx("./result/mlModelMetrics.xlsx", sheet = "rf_train")
trRF <- trRF[, c("drug", "geneFilterMethod", "pearsonCor")]
colnames(trRF) <- c("drug", "Type", "value")
pltTRRF <- create_plot(trRF)
## ------------------------------------------------------------
tsRF <- read_xlsx("./result/mlModelMetrics.xlsx", sheet = "rf_test")
tsRF <- tsRF[, c("drug", "geneFilterMethod", "pearsonCor")]
colnames(tsRF) <- c("drug", "Type", "value")
pltTSRF <- create_plot(tsRF)
## ------------------------------------------------------------
## ---- Elastic-Net training and test plot --------------------
trEN <- read_xlsx("./result/mlModelMetrics.xlsx", sheet = "glmnet_train")
trEN <- trEN[, c("drug", "geneFilterMethod", "pearsonCor")]
colnames(trEN) <- c("drug", "Type", "value")
pltTREN <- create_plot(trEN)
## ------------------------------------------------------------

tsEN <- read_xlsx("./result/mlModelMetrics.xlsx", sheet = "glmnet_test")
tsEN <- tsEN[, c("drug", "geneFilterMethod", "pearsonCor")]
colnames(tsEN) <- c("drug", "Type", "value")
pltTSEN <- create_plot(tsEN)
## ------------------------
pdf("result/az_Fig-2_ML_results.pdf", width = 8.3, height = 8.0)
print(ggarrange(pltTRRF, pltTSRF, pltTREN, pltTSEN,
  labels = c("A", "B", "C", "D"),
  ncol = 2, nrow = 2
))
dev.off()
print("done")
