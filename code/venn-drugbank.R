# library(ggpubr)
library(ggplot2)
library(ggvenn)
library(ggsci)
library(eulerr)
createVennDiagram <- function(l, title = NA) {
  p <- ggvenn(l, auto_scale = TRUE) +
    scale_fill_npg() +
    xlab("") +
    ylab("") +
    theme_void()
  return(p)
}


createEulerDiagram <- function(l, title = NA) {
  p <- plot(eulerr::euler(l),
    quantities = list(cex = .3),
    strips = list(cex = .3),
    labels = list(cex = .3),
    # lables = list(cex = .5, fontsize = 5),
    # fills = list(fill = c("#00A087", "#3C5488"), alpha = 0.85),
    fills = list(fill = c("#EA785D", "#88969C"), alpha = 1),
    main = title,
    fontsize = 10,
    edges = FALSE,
    cex.main = 0.25
  )
  print(p)
}

l <- list(
  # "Tozasertib" = c("RIPK1", "AURKA", "AURKB", "AURKC"),
  "Paclitaxel" = c(
    "BCL2", "TUBB1", "MAP4", "MAP2", "MAPT", "NR1I2",
    "CYP3A4", "CYP3A5", "CYP3A7", "CYP19A1", "CYP1B1",
    "CYP2C8", "ABCB11", "ABCB1", "ABCC1", "ABCC10", "SLCO1B3",
    "ABCC2"
  ),
  "Lapatinib" = c(
    "EGFR", "ERBB2",
    "CYP3A4", "CYP3A5", "CYP2C8", "CYP2C19",
    "ABCB1", "TAP1"
  ), # HER2",, "HER1","ERBB1",
  "Erlotinib" = c(
    "EGFR", "NR1I2", "CYP3A4", "CYP3A5", "CYP1A2", "CYP1A1",
    "CYP2D6", "CYP2C8", "CYP1B1", "UGT1A1", "ALB", "ORM1", "ABCG2", "ABCB1",
    "SLCO2B1"
  )
)
# maxn <- max(sapply(l, length))
# coln <- c("drug", paste0("gene_", seq_len(maxn)))
# out <- data.frame()
# l <- l[order(sapply(l, length))]
# for (i in names(l)) {
#   v <- c(l[[i]], rep("", maxn - length(l[[i]])))
#   out <- rbind(c(i, v), out)
# }
# colnames(out) <- coln
# write_xlsx(out, "./result/drugBank_gene_list.xlsx")

pdf("./result/venn-drugbank.pdf", height = 1, width = 3)
# pdf("./result/venn-drugbank.pdf", height = 3, width = 10)
# p <- list()
# vennList <- list()
# vennList2 <- list()
for (i in names(l)) {
  print(i)
  mat <- readRDS(sprintf("./data/drug_text-features/%s.rds", i))
  txtMineGenes <- unique(mat$Symbol)
  # p[[i]] <- createVennDiagram(vennList)
  # plot(createVennDiagram(vennList))
  vennList <- list("DrugBank" = unique(l[[i]]), "Text-mining genes" = txtMineGenes)
  vennList2 <- list("DrugBank" = unique(l[[i]]), "Text-mining genes" = txtMineGenes[seq_len(500)])
  createEulerDiagram(vennList, i)
  # createEulerDiagram(vennList2, i)
}
# createEulerDiagram(vennList[2], names(l)[2])
# createEulerDiagram(vennList2[2], names(l)[2])
# createEulerDiagram(vennList[3], names(l)[3])
# createEulerDiagram(vennList2[3], namthes(l)[3])
dev.off()
print("done")
