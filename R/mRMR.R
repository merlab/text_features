.fs.mRMR <- function(features, response, feature.count = 200, ...){
  args <- list(...)
  if ("mRMR.threadsCount" %in% names(args)){
    mRMRe::set.thread.count(args$mRMR.threadsCount)
  }
  dataMat <- mRMRe::mRMR.data(data=as.data.frame(cbind(response, features), stringAsFactor=FALSE))
  features <- mRMRe::mRMR.ensemble(data = dataMat, target_indices = 1, solution_count = 1, feature_count = feature.count)
  return(features@feature_names[unlist(features@filters)])
}