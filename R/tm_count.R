library("ggsci")
library("ggplot2")
library("gridExtra")
library(readr)
library(tidyr)
library(ggplot2)
library(Hmisc)
library(plyr)
library(RColorBrewer)
library(reshape2)
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# Plot number of text mining genes for each drug
counts <- readRDS("../data/tm_count.rds")
files <- list.files(path=genespath, full.names=FALSE, recursive=FALSE)
files <- file_path_sans_ext(files)
counts<-counts[counts$name %in% files,]
counts$count <- strtoi(counts$count)
counts <- counts[order(counts$count),]
counts$name <-factor(counts$name, levels = counts$name)
p <-ggplot(counts, aes(x=name, y=count, fill = "#4DBBD5FF")) + geom_bar(stat = "identity") + coord_flip()+ theme(legend.position="none")
p <- p + labs(title="Number of Text Mining Genes for each Drug",y="Number of Text Mining Genes", x = "Drug")
p
pdf("./figures/tm_count.pdf", width=6.5, height=9)

# Plot Pearson correlation
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

filename <- "paired_cor_train_CCLE_rf_class"
ds <- readRDS(sprintf("../data/dfs/%s.rds", filename))
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
dfcopy2$Type <- factor(as.character(dfcopy2$Type), levels = keycols2)

g <- ggplot(data = dfcopy2, aes(y = value, x = Type, fill = Type)) +
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .65) +
  geom_point(aes(y = value, color = Type), position = position_jitter(width = .15), size = .5, alpha = 0.8) +
  geom_boxplot(width = .1, guides = FALSE, outlier.shape = NA, alpha = 0.5) +
  expand_limits(x = 5.25) +
  guides(fill = FALSE) +
  guides(color = FALSE) +
  #scale_color_brewer(palette = "Spectral") +
  #scale_fill_brewer(palette = "Spectral") +
  #coord_flip() +
  theme_bw() +
  raincloud_theme + xlab("Feature Selection") + ylab("Pearson Correlation") + scale_color_npg()+ scale_fill_npg()
g
pdf(sprintf("./figures/rain/%s.pdf",filename), width=6.5, height=7)
print(g)
dev.off()

