library(ggplot2)
library(tools)

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