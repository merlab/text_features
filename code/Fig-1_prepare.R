source("./code/helper/summarizeData.R")
source("./code/helper/computeInteractionMatrix.R")
source("./code/helper/mRMR.R")
source("./code/helper/summarizeData_pharmacoGx.R")
library(PharmacoGx)
library(readxl)
library(writexl)
library(Xeva)
convertTissues <- function(v) {
  convL <- list(
    "Prostate" = c("Prostate", "Testis"),
    "Brain" = c("central_nervous_system", "autonomic_ganglia", "CNS", "brain", "Peripheral Nervous System"),
    "Kidney" = c("kidney", "Adrenal Gland"),
    "Bladder" = c("Bowel", "urinary_tract", "Bladder", "urine"),
    "Thyroid" = c("thyroid"),
    "Soft tissue" = c("soft_tissue", "Soft Tissue"),
    "Skin" = c("skin", "melanoma"),
    "Lung" = c("lung"),
    "Bone" = c("bone"),
    "Gynecological" = c("ovary", "pleura", "endometrium", "Pleura", "Cervix", "vulva", "vagina", "fallopian", "Uterus", "gestational"),
    "Pancreas" = c("pancreas"),
    "Breast" = c("breast"),
    "Liver" = c("liver"),
    "Digestive" = c("esophagus", "upper_aerodigestive_tract", "large_intestine", "biliary_tract", "stomach", "biliary"),
    # "Salivary Glands" = c("salivary_gland"),
    "Other" = c("haematopoietic_and_lymphoid_tissue", "Other", "Myeloid", "Lymphoid", "Head and Neck", "salivary_gland")
  )

  priotList <- c("Breast","Lung","Pancreas","Prostate","Gynecological",
                 "Brain","Kidney","Bladder","Thyroid",
  "Skin","Bone","Digestive","Liver",
  "Soft tissue",
  # "Salivary Glands",
  #"Head and Neck",
  "Other")
  # priotList <- rev(priotList)
  v <- tolower(v)
  for (i in names(convL)) {
    for (j in convL[[i]]) {
      v[grep(tolower(j), v, ignore.case = TRUE)] <- i
    }
  }
  v[v == "NA"] <- "Other"
  v[v == "na"] <- "Other"
  v <- factor(v, levels = priotList)
  return(v)
}
cleanNames <- function(x) {
  x <- na.omit(x)
  # x[is.na(x)] <- "NA"
  x[x == ""] <- "NA"
  return(x)
}
a <- readRDS("./data/ccle_ctrpv2_aac.rds")
a <- a[apply(a, 1, function(x) return(any(!is.na(x)))), ]
cclecells <- rownames(a)

a <- readRDS("./data/gdse_aac.rds")
a <- a[apply(a, 1, function(x) return(any(!is.na(x)))), ]
gdsecells <- rownames(a)

gdse <- readRDS("./data/PSet_GDSC2020.rds")
clinfo <- cellInfo(gdse)
# gdsecells <- rownames(clinfo)
gdsetissues <- clinfo$tissue
ccle <- readRDS("./data/CCLE-CTRPv2_Kallisto_0.46.1.rnaseq.rds")
clinfo <- cellInfo(ccle)
# cclecells <- rownames(clinfo)
ccletissues <- clinfo$ccle_primary_site



# cleaning
gdsecells <- unique(cleanNames(gdsecells))
gdsetissues <- cleanNames(gdsetissues)
gdsetissues <- convertTissues(gdsetissues)
cclecells <- unique(cleanNames(cclecells))
ccletissues <- cleanNames(ccletissues)
ccletissues <- convertTissues(ccletissues)
#
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
drugs <- list(
  "CCLE" = c(paste(rep("A", 38), seq_len(38)), paste(rep("B", 24), seq_len(24))),
  "GDSE" = c(paste(rep("A", 38), seq_len(38)))
)

cells <- list(
  "CCLE" = cclecells,
  "GDSE" = gdsecells
)


x <- as.data.frame(table(ccletissues))
y <- as.data.frame(table(gdsetissues))
colnames(x)[1] <- colnames(y)[1] <- "tissue"
colnames(x)[2] <- "n_cell_line_ccle"
colnames(y)[2] <- "n_cell_line_gdse"
#
m <- merge(x,y,by = "tissue")
write_xlsx(m, "./result/cell_line_per_tissue.xlsx")
saveRDS(drugs, "./drugs.rds")
saveRDS(cells, "./cells.rds")
saveRDS(ccletissues, "./ccletissues.rds")
saveRDS(gdsetissues, "./gdsetissues.rds")
