#' Get both molecular data and profile of a given drug from a PharmacoSet
#'
#' Given a PharmacoSet, this function will extract both molecular data
#' and drug response. This is done by running first summarizeMolecularProfiles
#' function and then summarizeSensitivityProfiles function. Note that this will
#' not work for perturbation type data.
#'
#' @examples
#' data(GDSCsmall)
#' df = summarizeData(pSet=GDSCsmall, mDataType = "rna", drug = "Sorafenib", sensitivity.measure="ic50_published")
#'
#' @param pSet \code{PharmacoSet} The PharmacoSet to summarize
#' @param mDataType \code{character} which one of the molecular data types
#' to use in the analysis, out of all the molecular data types available for the pset
#' for example: rna, rnaseq, snp
#' @param drug \code{character} The drug to be summarized. This will accept only one drug.
#' @param sensitivity.measure [character] which sensitivity sensitivity.measure to use? Use the
#'   sensitivityMeasures function to find out what measures are available for each PSet.
#' @param summary.stat.sensitivity \code{character} which summary method to use
#' if there are repeated cell line-drug experiments? Choices are "mean", "median", "first", or "last"
#' ... parameters for summarizeMolecularProfiles function
#'
#' @return An expressesion set object with the molecular data and drug profile.
#' 
#' @author Arvind Mer
#'
#' @export
summarizeData <- function(pSet, mDataType,
                 drug,
                 sensitivity.measure = NA,
                 summary.stat.sensitivity = c("mean", "median", "first", "last", "max", "min"),
                 ...)
{
  if (!drug %in% rownames(drugInfo(pSet)))
  {
    stop("drug not present, see drugInfo(pSet)")
  }
  
  if (!is.na(sensitivity.measure))
  {
    if (sensitivity.measure %in% sensitivityMeasures(pSet))
    {
      asm <- sensitivity.measure
    } else
    {
      stop("sensitivity.measure not present in pSet, see sensitivityMeasures(pSet)")
    }
  } else {
    asm <- sensitivityMeasures(pSet)
  }
  
  exp <- summarizeMolecularProfiles(pSet = pSet, mDataType = mDataType, ...)
  clid <- exp$cellid
  for (s in asm)
  {
    sm <- summarizeSensitivityProfiles(
      pSet = pSet,
      sensitivity.measure = s,
      cell.lines = clid,
      drugs = drug,
      summary.stat = summary.stat.sensitivity[1],
      fill.missing = TRUE,
      verbose = TRUE
    )
    
    if(is(exp, "SummarizedExperiment"))
    {
      colData(exp)[, s] <- sm[drug, clid]
    } else { pData(exp)[, s] <- sm[drug, clid]} 
  }
  
  return(exp)
}
