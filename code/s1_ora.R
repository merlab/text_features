library(fgsea)
library(ggplot2)
library(writexl)
filterPathMat <- function(x) {
  x <- as.data.frame(x)
  x <- x[x$pval <= 0.05, ]
  x <- x[x$padj <= 0.05, ]
  x$overlapGenes <- NULL
  return(x)
}
# load pathways
c2 <- gmtPathways("./data/pathways/c2.all.v2023.2.Hs.symbols.gmt")
c2allGenes <- unique(unlist(c2))
go <- gmtPathways("./data/pathways/c5.go.v2023.2.Hs.symbols.gmt")
goAllGenes <- unique(unlist(go))
# find text mining gnges
drugs <- list.files("./data/drug_text-features/", pattern = "*.rds")
drugs <- gsub("\\.rds", "", drugs)
#### ADD text mining genes to the universes for the pathways
for (drug in drugs) {
  genes <- readRDS(sprintf("./data/drug_text-features/%s.rds", drug))$Symbol
  c2allGenes <- unique(c(c2allGenes, genes))
  goAllGenes <- unique(c(goAllGenes, genes))
}
#### fora ORA analysis
c2res <- list()
gores <- list()
for (drug in drugs) {
  genes <- readRDS(sprintf("./data/drug_text-features/%s.rds", drug))$Symbol
  c2res[[drug]] <- fora(c2, genes = genes, universe = c2allGenes, minSize = 20)
  gores[[drug]] <- fora(go, genes = genes, universe = goAllGenes, minSize = 20)
}
gores <- lapply(gores, filterPathMat)
c2res <- lapply(c2res, filterPathMat)
write_xlsx(gores, "./result/textMining_ORA_GO.xlsx")
write_xlsx(c2res, "./result/textMining_ORA_C2.xlsx")

print("done")
