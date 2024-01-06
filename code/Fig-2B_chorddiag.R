# loading libraries
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
library(tidyr)
library(readxl)
library(awtools)

##################
### parameters ###
##################
# if gene is in the top x of the drug connection is shown
n_top_drug <- 5
# x genes with most connections shows
n_genes_overall <- 15
# x drugs with most connections shown
n_drugs_overall <- 38

# limit to drugs that are approved based on the
lim_drug_bank <- TRUE


#################################
### preparing plotting dataset###
#################################
# gene and drug sorting function in the chord diag
#   genes/drugs with the most number of connections are at the center
#   genes/drug with least number of connections are at the periphery
new_sorting <- function(v) {
  o <- rep(NA, length(v))
  mid <- length(v) %/% 2 + 1
  for (i in 1:length(v)) {
    n <- v[i]
    if (i == 1) {
      o[mid] <- n
    } else if (i %% 2 == 0) {
      o[mid - (i %/% 2)] <- n
    } else {
      o[mid + (i %/% 2)] <- n
    }
  }
  return(o)
}

# loading the text mining data
dir <- "./data/drug_text-features"
df <- data.frame()
for (i in list.files(dir, patter = "*.rds")) {
  drug <- gsub("\\.rds", "", i)
  t <- readRDS(sprintf("%s/%s", dir, i))
  t <- t[t$FDR < 0.05, ]
  t <- t[order(t$FDR, decreasing = FALSE), ]
  t <- t[seq_len(min(n_top_drug, nrow(t))), ]
  t <- data.frame(from = drug, to = t$Symbol, w = 1)
  df <- rbind(t, df)
}
# make matrix of drugs vs genes
mat <- pivot_wider(df, names_from = to, values_from = w, values_fill = 0)
mat <- as.data.frame(mat)
# fix the rownames
rownames(mat) <- mat$from
mat$from <- NULL

# filter the drugs by the drugs in the approved stat
if (lim_drug_bank == TRUE) {
  filter <- read_xlsx("~/text_features_all/text_features/result/Supplementary_Data_File-1_txtMiningGeneList.xlsx")
  filter <- filter$drug[filter$approved == TRUE]
  mat <- mat[rownames(mat) %in% filter, ]
}

# sort for top genes
mat <- mat[order(rowSums(mat), decreasing = TRUE), ]
# sort for top drugs
mat <- mat[, order(colSums(mat), decreasing = TRUE), ]

# choose the top x genes and top y drugs as degined in parameters
mat <- mat[1:n_drugs_overall, 1:n_genes_overall]

# sort for top genes again
mat <- mat[order(rowSums(mat), decreasing = TRUE), ]
# sort for top drugs again
mat <- mat[, order(colSums(mat), decreasing = TRUE), ]

mat <- as.matrix(mat)
mat <- na.omit(mat)
# for coloring and direction of chord diag
mat <- t((mat))

# the gene and drug sorting
mat <- mat[new_sorting(rownames(mat)), ]
mat <- mat[, new_sorting(colnames(mat))]

# color formatting for the drugs
drug_cols <- unique(c(
  "#ff420e", "#ffd320", "#7e0021", "#83caff", "#314004", "#aecf00", "#4b1f6f",
  "#ff950e", "#c5000b", "#3366cc", "#dc3912", "#109618", "#990099", "#0099c6",
  "#dd4477", "#66aa00", "#b82e2e", "#316395", "#2A363B", "#019875", "#FECEA8",
  "#FF847C", "#E84A5F", "#C0392B", "#96281B", "#F7DC05", "#EC0B88", "#5e35b1",
  "#f9791e", "#3dd378", "#c6c6c6", "#444444", "#017a4a", "#FFCE4E", "#3d98d3",
  "#ff363c", "#7559a2", "#794924", "#8cdb5e", "#d6d6d6", "#fb8c00"
))[1:ncol(mat)]



# color used for the genes
gene_cols <- rep("#57635F", nrow(mat))

# all colors used in the plot
cols <- c(drug_cols, gene_cols)

# open pdf
pdf("./result/Fig-2B_chorddiag.pdf", width = 8, height = 8)

par(mar = c(0.5, 3, 1, 0), xpd = NA)

# orientation parameters for the plot to be horizontal
circos.par(
  start.degree = -90
  # space between the plot
  , track.margin = c(0.001, 0.0001)
)
chordDiagramFromMatrix(
  # the input data
  t(mat),
  keep.diagonal = TRUE,
  # required parameters from the code i used to correct name directions
  annotationTrack = "grid", preAllocateTracks = 1,
  ### beautification parameters ###
  # for colors of the items on the diagram
  grid.col = cols,
  # for colors of the link
  row.col = drug_cols,
  # line border type
  link.lty = 1,
  # line border thickness
  link.lwd = 0.002,
  # it has to be 1 for the borders to be drawn
  link.border = 0,
  # transparency of the lanes
  transparency = 0.5,
)

# I have no clue about how it works
# it throws out many errors but it corrects the names on the chord diag
# https://bioinfo4all.wordpress.com/2021/03/13/tutorial-7-how-to-do-chord-diagram-using-r/
# track index of 2 is needed to get the names close to the circle
circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  sector.name <- get.cell.meta.data("sector.index")
  # number after ylim determines how far the labels are from the plot
  circos.text(mean(xlim), ylim[1] + 1.5, sector.name,
    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
    # relative font size of the names
    cex = 1.35
  )
}, bg.border = NA)

circos.clear()
# dev.off()
print("done")
