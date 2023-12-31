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
df <- df[df$drugName %in% drugs, ]

# common_drugs <- intersect(f, drugs)
# stop()
min(df$n_common_cells[df$drugName %in% f])
min(df$n_cells_CCLE[df$drugName %in% f])
min(df$n_cells_GDSE[df$drugName %in% f])
# table(df$n_cells_CCLE > 100) # & df$n_cells_GDSE > 100)
# table(df$n_cells_GDSE > 250)
# table(df$n_common_cells > 250)
fda <- read.csv("./data/FDA-approved-drug-list.csv")
drugNames <- tolower(unique(c(fda$Generic.Name, fda$Brand.Name)))

#
df$isFDAapproved <- sapply(df$drugName, function(x) {
  return(any(grepl(x, drugNames)))
})
table(df$isFDAapproved & df$n_cells_GDSE > 0)
table(df$isFDAapproved & df$n_cells_CCLE > 0)
# grep(df$drugName[1], drugNames)
# table(tolower(df$drugName) %in% drugNames & df$n_cells_CCLE > 0)
# table(tolower(df$drugName) %in% drugNames & df$n_cells_GDSE > 0)
