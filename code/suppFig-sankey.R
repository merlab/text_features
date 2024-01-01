library(tidyr)
library(readxl)
library(awtools)
library(ggplot2)
library(ggsankey)
library(ggalluvial)
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


##################
### parameters ###
##################
# limit to drugs that are approved based on the
lim_drug_bank <- TRUE
n_top_drug <- 5
#
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


mat <- pivot_wider(df, names_from = to, values_from = w, values_fill = 0)
mat <- as.data.frame(mat)
# fix the rownames
rownames(mat) <- mat$from
mat$from <- NULL

# drugs <- new_sorting(rownames(mat))
# genes <- new_sorting(colnames(mat))
#####
v <- setNames(colSums(mat), colnames(mat))
v <- v[v > 0]
df <- df[df$to %in% names(v), ]
mat <- mat[, colnames(mat) %in% names(v)]
# #####
# v <- setNames(rowSums(mat), rownames(mat))
# v <- v[v > 2]
# df <- df[df$from %in% names(v), ]
# mat <- mat[rownames(mat) %in% names(v), ]
# #####
#####
#####
v <- setNames(colSums(mat), colnames(mat))
genes <- names(sort(v, decreasing = TRUE))
#####
v <- setNames(rowSums(mat), rownames(mat))
drugs <- names(sort(v, decreasing = TRUE))
#####
#####
#####
#####
# color formatting for the drugs
drug_cols <- unique(c(
  "#ff420e", "#ffd320", "#7e0021", "#83caff", "#314004", "#aecf00", "#4b1f6f",
  "#ff950e", "#c5000b", "#3366cc", "#dc3912", "#109618", "#990099", "#0099c6",
  "#dd4477", "#66aa00", "#b82e2e", "#316395", "#2A363B", "#019875", "#FECEA8",
  "#FF847C", "#E84A5F", "#C0392B", "#96281B", "#F7DC05", "#EC0B88", "#5e35b1",
  "#f9791e", "#3dd378", "#c6c6c6", "#444444", "#017a4a", "#FFCE4E", "#3d98d3",
  "#ff363c", "#7559a2", "#794924", "#8cdb5e", "#d6d6d6", "#fb8c00"
))[seq_along(unique(df$from))]
gene_cols <- rep(c("#889296", "#D5D5D5"), length(unique(df$to)) %/% 2)
# if number of genes is odd add a filler
if (length(unique(df$to)) %% 2 > 0) gene_cols <- c(gene_cols, "#889296")
cols <- c(drug_cols, gene_cols)
#####
#####
### formatting for geom sankey
df$from <- factor(df$from, levels = rev(drugs))
df$to <- factor(df$to, levels = rev(genes))
df$x <- 1
df$next_x <- 2
df$node <- as.numeric(df$from)
df$next_node <- as.numeric(df$to)
df2 <- data.frame(from = df$to, to = NA, w = 1, x = 2, next_x = 3, node = as.numeric(df$to), next_node = NA)
df <- rbind(df, df2)
#####
#####
p <- ggplot(
  df,
  aes(
    x = x, next_x = next_x, node = node, next_node = next_node,
    label = from, fill = from
  )
) +
  geom_sankey(aes(node.fill = from), alpha = 0.75) +
  geom_sankey_label(color = 1, fill = "white") +
  scale_fill_manual(values = cols) +
  theme_void() +
  theme(legend.position = "none")

# p <- ggplot(df, aes(axis1 = from, axis2 = to, y = w)) +
#   geom_alluvium(aes(fill = from)) +
#   geom_stratum(aes(fill = from)) +
#   geom_text(
#     stat = "stratum",
#     aes(label = after_stat(stratum))
#   ) +
#   scale_fill_manual(values = cols) +
#   # scale_fill_manual(values = drug_cols) +
#   theme_void() +
#   theme(legend.position = "none")
pdf("./result/sankey.pdf", height = 10.5 * 3.5, width = 7.5 * 2)
plot(p)
dev.off()
print("done")
