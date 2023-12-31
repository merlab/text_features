library(tidyxl)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(viridis)
source("./code/helper/pathwayName_cleaner.R")
dotplotORA <- function(m, title = NA) {
  m$size <- as.numeric(m$size)
  # <- m[m$size < 1250, ]
  m$overlap <- as.numeric(m$overlap)
  m$padj <- as.numeric(m$padj)
  m$pval <- as.numeric(m$pval)
  m$pathway <- cleanNames(m$pathway)
  #
  # torm <- rmExtra(m$pathway)
  # if (length(torm) > 0) {
  #   m <- m[-torm, ]
  # }
  #
  m <- rmExtra(m)
  chars <- unlist(sapply(m$pathway, nchar))
  m <- m[chars < 60, ]
  #
  # m <- m[m$padj <= 0.15, ]
  m$percentage <- m$overlap / m$size * 100
  # m <- m[m$pval <= 0.05, ]
  m$logFDR <- -log10(m$padj)
  m$logPVL <- -log10(m$pval)
  m$logFDR[m$logFDR > 300] <- 300
  m$logFDR[m$padj == 0] <- 300
  # print(apply(m, 2 , summary))
  # m <- m[seq_len(min(nrow(m), 15)), ]
  m <- m[order(m$percentage, decreasing = FALSE), ]
  if (nrow(m) == 0) {
    return(NULL)
  } else {
    m <- m[!duplicated(m$pathway), ]
    m$pathway <- factor(m$pathway, levels = m$pathway)
    p <- ggplot(m, aes(x = percentage, y = pathway, size = size, color = logPVL)) + # logFDR
      geom_point() +
      scale_color_viridis() +
      scale_size_continuous(breaks = seq(0, 1500, 250), range = c(2, 5)) +
      xlab("Percentage Overlap (%)") +
      ylab("") +
      labs(colour = expression("-log"[10] * "FDR"), size = "Pathway Size") +
      # scale_x_continuous(breaks = seq(0, 100, 25), limits = c(0, 100)) +
      ggtitle(title) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 10),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent", color = NA) # get rid of legend panel bg
      )
    return(p)
  }
}

formats <- xlsx_formats("./result/textMining_ORA_C2.xlsx")
formats$local$fill$patternFill$fgColor$rgb
cells <- xlsx_cells("./result/textMining_ORA_C2.xlsx")
sheets <- unique(cells$sheet)
subcellsC2 <- cells[
  which(cells$local_format_id %in% which(!is.na(formats$local$fill$patternFill$fgColor$rgb))) - 1,
  c("sheet", "row")
]

formats <- xlsx_formats("./result/textMining_ORA_GO.xlsx")
formats$local$fill$patternFill$fgColor$rgb
cells <- xlsx_cells("./result/textMining_ORA_GO.xlsx")
sheets <- unique(cells$sheet)
subcellsGO <- cells[
  which(cells$local_format_id %in% which(!is.na(formats$local$fill$patternFill$fgColor$rgb))) - 1,
  c("sheet", "row")
]

p <- list()
for (i in sheets) {
  c2Mat <- as.data.frame(readxl::read_xlsx("./result/textMining_ORA_C2.xlsx", sheet = i))
  if (!any(subcellsC2$sheet %in% i)) next()
  print(i)

  rows <- unique(subcellsC2$row[subcellsC2$sheet %in% i])
  c2SubMat <- c2Mat[rows, ]
  # c2Mat <- c2Mat[-rows, ]
  # c2SubMat <- rbind(c2Mat, c2Mat[seq_len(15 - nrow(c2Mat)), ])


  goMat <- as.data.frame(readxl::read_xlsx("./result/textMining_ORA_GO.xlsx", sheet = i))
  if (any(subcellsGO$sheet %in% i)) {
    print(i)
    rows <- unique(subcellsGO$row[subcellsGO$sheet %in% i])
    goSubMat <- goMat[rows, ]
    oraMat <- rbind(c2SubMat, goSubMat)
  } else {
    oraMat <- c2SubMat
  }
  if (i == "Paclitaxel") {
    oraMat <- rbind(oraMat, c2Mat[grep("tubu", c2Mat$pathway, ignore.case = TRUE), ])
  }
  # goMat <- goMat[-rows, ]
  # goSubMat <- rbind(goMat, goMat[seq_len(15 - nrow(goMat)), ])
  p[[i]] <- dotplotORA(oraMat, i)
}
pdf("./result/ORA_res.pdf", height = 10, width = 12)
plot(ggarrange(
  plotlist = p[1:2], nrow = 2, ncol = 1,
  labels = paste0("    (", LETTERS[1:2], ")"),
  common.legend = TRUE, align = "v", legend = "right"
))
plot(ggarrange(
  plotlist = p[3:4], nrow = 2, ncol = 1,
  labels = paste0("    (", LETTERS[3:4], ")"),
  common.legend = TRUE, align = "v", legend = "right"
))
dev.off()
