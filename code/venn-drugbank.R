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

l <- list(
  "Tozasertib" = c("RIPK1", "AURKA", "AURKB", "AURKC"),
  "Paclitaxel" = c("BCL2", "TUBB1", "MAP4", "MAP2", "MAPT", "NR1I2"),
  "Lapatinib" = c("EGFR", "HER2", "ERBB2", "HER1", "ERBB1"), #
  "Erlotinib" = c("EGFR", "NR1I2")
)

p <- list()
for (i in names(l)) {
  mat <- readRDS(sprintf("./data/drug_text-features/%s.rds", i))
  txtMineGenes <- mat$Symbol
  vennList <- list("DrugBank" = l[[i]], "Genie" = txtMineGenes)
  p[[i]] <- createVennDiagram(vennList)
}

pdf("./result/venn-drugbank.pdf", height = 6, width = 6)
plot(ggarrange(
  plotlist = p, nrow = 2, ncol = 2,
  labels = paste0("(", LETTERS[1:4], ") ", names(l))
))
dev.off()
