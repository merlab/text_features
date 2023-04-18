library(ggplot2)
library(ggpubr)
#library(gridExtra)
#library(readr); library(tidyr)
#library(tools); library(Hmisc) ; library(plyr); library(RColorBrewer)
#library(reshape2); ; library(ggrepel)

#source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
source("R/geom_flat_violin.R")

raincloud_theme = theme(
  text = element_text(size = 10),
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.title=element_text(size=16),
  legend.text=element_text(size=16),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 16),
  panel.border = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

##-------------------------------------------------------------

create_plot <- function(df, title=NA)
{
  fetMeth <- colnames(df)[colnames(df)!="drug"]

  cl <- c("#B09C85FF", "#91D1C2FF", "#4DBBD5FF","#3C5488FF", "#00A087FF","#DC0000FF")
  names(cl) <- fetMeth

  my_comparisons <- lapply(rev(fetMeth[1:(length(fetMeth)-1)]),
                           function(i) c(i, fetMeth[length(fetMeth)]))

  df$drug <- rownames(df)

  dft <- reshape2::melt(df, id.vars='drug', variable.name = "Type")

  g <- ggplot(data = dft, aes(y = value, x = Type, fill = Type)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .95) +
    geom_point(aes(y = value, color = Type),
               position = position_jitter(width = .15), size = .5, alpha = 0.45) +
    geom_boxplot(width = .1, #guides = "none",
                 outlier.shape = NA, alpha = 0.5) +
    expand_limits(x = 5.25) + guides(fill = "none") + guides(color = "none") +
    theme_bw() +theme(plot.title = element_text(hjust = 0.5))+
    raincloud_theme + xlab("")+ #xlab("Feature Selection Method") +
    ylab("Pearson correlation")

  g <- g + scale_color_manual(values = cl) + scale_fill_manual(values = cl)
  g <- g + stat_compare_means(comparisons = my_comparisons, method = "t.test", paired=TRUE,
                              #label.y = 0.7,
                              size = 2.5)
  g <- g + scale_y_continuous(breaks =seq(0,1,0.1))
  g + theme(axis.text.x = element_text(face="bold", #color="#993333",
                                      size=8, angle=-0),
           axis.text.y = element_text(face="bold", #color="#993333",
                                      size=9, angle=0),
           axis.title.y = element_text(face="bold", size=10))
}

##--------------------------------------------------------------
##---- Random-forest training and test plot --------------------
trRF <- readRDS('result/train_test/train_RandomForest_Reg.rds')
pltTRRF <- create_plot(trRF)
print(pltTRRF+ggtitle("Random Forest training"))

tsRF <- readRDS('result/train_test/test_RandomForest_Reg.rds')
pltTSRF <- create_plot(tsRF)
print(pltTSRF+ggtitle("Random Forest test"))

##------------------------------------------------------------
##---- Elastic-Net training and test plot --------------------
trEN <- readRDS('result/train_test/train_ElasticNet_Reg.rds')
pltTREN <- create_plot(trEN)
print(pltTREN+ggtitle("Elastic-Net training"))


tsEN <- readRDS('result/train_test/test_ElasticNet_Reg.rds')
pltTSEN <- create_plot(tsEN)
print(pltTSEN+ggtitle("Elastic-Net test"))

##------------------------
pdf("result/Fig-2_ML_results.pdf", width = 8.3, height = 8.0)
print(ggarrange(pltTRRF, pltTSRF, pltTREN, pltTSEN,
          #labels = c("A", "B", "C"),
          ncol = 2, nrow = 2))
dev.off()



