library(ggpubr)
library(ggplot2)
library(ggvenn)
library(ggsci)
createVennDiagram <- function(l, title = NA) {
  p <- ggvenn(l) +
    scale_fill_npg() +
    xlab("") +
    ylab("") +
    theme_void()
  return(p)
}
createPieChart <- function(v, title = NA) {
  df <- as.data.frame(table(v))
  colnames(df) <- c("tissue", "freq")
  df <- df[order(df$freq, decreasing = TRUE), ]
  df$tissue <- factor(df$tissue, levels = df$tissue)
  p <- ggplot(df, aes(x = "", y = freq, fill = tissue)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    scale_fill_manual(
      name = "",
      values = c(
        "grey75",
        "#FED439FF", "#709AE1FF", "magenta", "#D2AF81FF", "#FD7446FF",
        "#D5E4A2FF", "#197EC0FF", "#F05C3BFF", "#46732EFF", "#71D0F5FF",
        "#370335FF", "#075149FF", "#C80813FF", "#91331FFF", "#1A9993FF",
        "#FD8CC1FF"
        # "#8A9197FF",
      )
    ) +
    coord_polar("y", start = 0) +
    xlab("") +
    ylab("") +
    theme_void()
  return(p)
}

cells <- readRDS("./cells.rds")
drugs <- readRDS("./drugs.rds")
ccletissues <- readRDS("./ccletissues.rds")
gdsetissues <- readRDS("./gdsetissues.rds")
pdf("./result/fig-1_venn.pdf", height = 5, width = 10)
venn1 <- createVennDiagram(cells)
venn2 <- createVennDiagram(drugs)
plot(ggarrange(venn1, venn2, nrow = 1, ncol = 2, common.legend = TRUE))
pie1 <- (createPieChart(ccletissues, "CCLE"))
pie2 <- (createPieChart(gdsetissues, "GDSE"))
plot(ggarrange(pie1, pie2, labels = c("CCLE", "GDSE"), nrow = 1, ncol = 2, common.legend = TRUE, legend = "right"))
dev.off()
print("done")
