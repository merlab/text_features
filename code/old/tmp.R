library(readxl)
removeCombTherapies <- function(v) {
  v <- v[!grepl(":", v)]
  v <- v[!grepl("^ML", v)]
  v <- v[!grepl("^MK", v)]
  v <- v[!grepl("^VER", v)]
  v <- v[!grepl("^VU", v)]
  v <- v[!grepl("^SKI", v)]
  v <- v[!grepl("^SJ", v)]
  v <- v[!grepl("[0-9]", v)]
  v <- v[!grepl("^N-", v)]
  v <- v[!grepl("^\\(-\\)-", v)]
  v <- v[!grepl("-leu-", v)]
  return(v)
}
df <- read_xlsx("./result/common_cell_lines_per_drug.xlsx")
df$drugName <- tolower(df$drugName)
df <- as.data.frame(df)
f <- list.files("./data/drug_text-features/")
f <- tolower(gsub("\\.rds", "", f))


drugs <- tolower(removeCombTherapies(df$drugName))

common_drugs <- intersect(f, drugs)
# stop()
df <- df[df$drugName %in% drugs, ]
min(df$n_common_cells[df$drugName %in% f])
min(df$n_cells_CCLE[df$drugName %in% f])
min(df$n_cells_GDSE[df$drugName %in% f])
table(df$n_cells_CCLE > 100) # & df$n_cells_GDSE > 100)
# table(df$n_cells_GDSE > 250)
# table(df$n_common_cells > 250)
