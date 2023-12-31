library(writexl)
library(readxl)
df <- read_xlsx("./result/common_cell_lines_per_drug.xlsx")

f <- list.files("./data/drug_text-features/")
f <- gsub("\\.rds", "", f)
write_xlsx(df[df$drugName %in% f, ], "./tmp.xlsx")
