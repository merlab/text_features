source("./code/helper/summarizeData.R")
source("./code/helper/computeInteractionMatrix.R")
source("./code/helper/mRMR.R")
source("./code/helper/summarizeData_pharmacoGx.R")
library(PharmacoGx)
library(readxl)
library(Xeva)
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
createPieChart <- function(v, title = NA) {
  df <- as.data.frame(table(v))
  colnames(df) <- c("tissue", "freq")
  df <- df[order(df$freq, decreasing = TRUE), ]
  df$tissue <- factor(df$tissue, levels = df$tissue)
  p <- ggplot(df, aes(x = "", y = freq, fill = tissue)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    scale_fill_simpsons(name = "") +
    coord_polar("y", start = 0) +
    xlab("") +
    ylab("") +
    theme_void()
  return(p)
}
createBarChartDf <- function(v) {
}
convertTissues <- function(v) {
  convL <- list(
    "Prostate" = c("Prostate", "Testis"),
    "Brain" = c("central_nervous_system", "autonomic_ganglia", "CNS", "brain", "Peripheral Nervous System"),
    "Kidney" = c("kidney", "Adrenal Gland"),
    "Bladder" = c("Bowel", "urinary_tract", "Bladder", "urine"),
    "Thyroid" = c("thyroid"),
    "Soft_tissue" = c("soft_tissue", "Soft Tissue"),
    "Salivary Glands" = c("salivary_gland"),
    "Head and Neck" = c("Head and Neck"),
    "Skin" = c("skin"),
    "Lung" = c("lung"),
    "Bone" = c("bone"),
    "Gynecological" = c("ovary", "pleura", "endometrium", "Pleura", "Cervix", "vulva", "vagina", "fallopian", "Uterus", "gestational"),
    "Pancreas" = c("pancreas"),
    "Breast" = c("breast"),
    "Liver" = c("liver"),
    "Digestive" = c("esophagus", "upper_aerodigestive_tract", "large_intestine", "biliary_tract", "stomach", "biliary"),
    "Other" = c("haematopoietic_and_lymphoid_tissue", "Other", "Myeloid", "Lymphoid")
  )
  v <- tolower(v)
  for (i in names(convL)) {
    for (j in convL[[i]]) {
      v[grep(tolower(j), v, ignore.case = TRUE)] <- i
    }
  }
  v[v == "NA"] <- "Other"
  v[v == "na"] <- "Other"
  return(v)
}
cleanNames <- function(x) {
  x <- na.omit(x)
  # x[is.na(x)] <- "NA"
  x[x == ""] <- "NA"
  return(x)
}
gdse <- readRDS("./data/PSet_GDSC2020.rds")
clinfo <- cellInfo(gdse)
gdsecells <- rownames(clinfo)
gdsetissues <- clinfo$tissue
ccle <- readRDS("./data/CCLE-CTRPv2_Kallisto_0.46.1.rnaseq.rds")
clinfo <- cellInfo(ccle)
cclecells <- rownames(clinfo)
ccletissues <- clinfo$ccle_primary_site



# cleaning
gdsecells <- cleanNames(gdsecells)
gdsetissues <- cleanNames(gdsetissues)
gdsetissues <- convertTissues(gdsetissues)
cclecells <- cleanNames(cclecells)
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
df <- read_xlsx("./result/common_cell_lines_per_drug.xlsx")
ccledrugs <- df$drugName[df$n_cells_CCLE > 100]
gdsedrugs <- df$drugName[df$n_cells_GDSE > 0]

ccledrugs <- tolower(removeCombTherapies(ccledrugs))
gdsedrugs <- tolower(removeCombTherapies(gdsedrugs))



txtMineDrugs <- list.files("./data/drug_text-features/")
txtMineDrugs <- tolower(gsub("\\.rds", "", txtMineDrugs))
common_drugs <- intersect(ccledrugs, gdsedrugs)
excluded_drugs <- common_drugs[!common_drugs %in% txtMineDrugs]
txtMineDrugs[!txtMineDrugs %in% common_drugs]
ccledrugs <- ccledrugs[!ccledrugs %in% excluded_drugs]
gdsedrugs <- gdsedrugs[!gdsedrugs %in% excluded_drugs]


drugs <- list("GDSE" = gdsedrugs, "CCLE" = ccledrugs)
cells <- list("GDSE" = gdsecells, "CCLE" = cclecells)
pdf("./result/fig-1_venn.pdf", height = 5, width = 10)
venn1 <- createVennDiagram(cells)
venn2 <- createVennDiagram(drugs)
plot(ggarrange(venn1, venn2, nrow = 1, ncol = 2, common.legend = TRUE))
pie1 <- (createPieChart(ccletissues, "CCLE"))
pie2 <- (createPieChart(gdsetissues, "GDSE"))
plot(ggarrange(pie1, pie2, labels = c("CCLE", "GDSE"), nrow = 1, ncol = 2, common.legend = TRUE, legend = "right"))
dev.off()
print("done")
