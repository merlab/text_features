library("ggsci")
library("ggplot2")
library("gridExtra")
library(readr)
library(tidyr)
library(ggplot2)
library(tools)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(ggpubr)
library(ggrepel)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# Plot Pearson correlation

# Read data
filename <- "paired_cor_train_CCLE_rf"
ds <- readRDS(sprintf("../data/dfs/%s.rds", filename))

# Prepare data for plotting
keycols <- c("var-100", "var-500", "L1000-tm", "L1000","cor-500", "text-mining-500")
keycols2 <- c("var-100", "var-500", "L1000-tm", "L1000","cor-500", "text-mining")
df <- data.frame(matrix(unlist(ds), ncol = max(lengths(ds)), byrow = TRUE))
names(df) <- names(ds[[which(lengths(ds)>0)[1]]])
df <- t(df)
colnames(df) <-names(ds)
dfcopy <- as.data.frame(df)[,keycols]
dfcopy$drug <- rownames(dfcopy)
colnames(dfcopy)[6] <- "text-mining"
dfcopy2 <- reshape2::melt(dfcopy, id.vars='drug', variable.name = "Type")
dfcopy2 <- dfcopy2[dfcopy2$drug %in% least$name,]
dfcopy2$Type <- factor(as.character(dfcopy2$Type), levels = keycols2)

my_comparisons <- list(c("cor-500", "text-mining"),c("L1000", "text-mining"),c("L1000-tm", "text-mining"), c("var-500", "text-mining"), c("var-100", "text-mining") )

# Raincloud Plot
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
g <- ggplot(data = dfcopy2, aes(y = value, x = Type, fill = Type)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .65) +
  geom_point(aes(y = value, color = Type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  theme_bw() +
  raincloud_theme + xlab("Feature Selection Method") + ylab("Pearson Correlation") + scale_color_npg()+ scale_fill_npg() +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", paired=TRUE)
g
pdf(sprintf("./figures/rain/%s.pdf",filename), width=6.5, height=7)
print(g)
dev.off()

# Line plot
plt <- ggplot(dfcopy2, aes_string(x="Type", y="value", fill = "Type")) +
  geom_line(aes(group=drug), color="#1f78b4", alpha = 0.35)
plt <- plt + theme(legend.position="none", axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                   axis.text = element_text(size = 14),panel.background = element_blank(),axis.line.x = element_line(color="black"),
                   axis.line.y = element_line(color="black")) + 
  xlab("Feature Selection Method") + ylab("Pearson Correlation") 
plt <- plt +     geom_text_repel(
  data = dfcopy2[dfcopy2$Type == "text-mining",],
  aes(label = drug),
  size = 2,
  direction = "y",
  xlim = c(6.3, NA),
  hjust = 0,
  force = 1,
  box.padding = .4,
  segment.size = 0.3,
  segment.color = "grey"
) + scale_x_discrete(
  expand = expansion(mult = c(.1, .2))
) + geom_point(color="#2166ac")
plt
pdf(sprintf("./figures/line/%s.pdf",filename), width=6.5, height=5)
print(plt)
dev.off()