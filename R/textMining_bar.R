library(ggplot2)
# loading the text mining data
dir <- "./data/drug_text-features"
plotdf <- data.frame()
for (i in list.files(dir, patter = "*.rds")) {
    drug <- gsub("\\.rds", "", i)
    x <- readRDS(sprintf("%s/%s", dir, i))
    x <- x[x$FDR < 0.05, ]
    x <- x[order(x$FDR, decreasing = FALSE), ]
    plotdf <- rbind(c(drug = drug, nGene = nrow(x)), plotdf)
}
colnames(plotdf) <- c("drug", "nGene")
plotdf$nGene <- as.numeric(plotdf$nGene)
p <- ggplot(plotdf, aes(x = reorder(drug, nGene), y = nGene)) +
    geom_bar(stat = "identity", color = NA, width = .8, fill = "#A4BC92") +
    # "grey50"   "#BA90C6" "#F4B183" "#804674" "#A86464"
    coord_flip() +
    scale_y_continuous(expand = c(0, 0)) +
    xlab("") +
    ylab("# of Genes") +
    theme_classic() +
    theme(
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
    )
pdf("./result/txtMineBar.pdf", height = 8, width = 4)
plot(p)
dev.off()
print("done")
