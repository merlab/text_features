library(Xeva); library(Biobase)
fl <- "~/CXP/XG/PDXE_xevaSet/data/PDXE_xevaSet/PDXE_XevaSet_All.rds"
df <- readRDS(fl)
##-----------------------------
##-------RNAseq ---------------
md <- df@molecularProfiles$RNASeq
saveRDS(md, "data/PDXE/PDXE_RNAseq.rds")
for(i in unique(md$tissue))
{
  mdt <- md[, md$tissue==i]
  if(dim(mdt)[2] >20 )
  {
    saveRDS(mdt, sprintf("data/PDXE/PDXE_RNAseq_%s.rds",i))
  }
}

##---------------------------------
##-------microarray ---------------
md <- df@molecularProfiles$microArray
saveRDS(md, "data/PDXE/PDXE_microArray.rds")
for(i in unique(md$tumor.type))
{
  mdt <- md[, md$tumor.type==i]
  if(dim(mdt)[2] >40 )
  {
    saveRDS(mdt, sprintf("data/PDXE/PDXE_microArray_%s.rds",i))
  }
}
