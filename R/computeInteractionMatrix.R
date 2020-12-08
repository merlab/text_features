
distancePointLine <-
  function(x, #x-coordinate of point
           y, #y-coordinate of point
           a, #coefficient in line equation ax + by + c = 0
           b, #coefficient in line equation ax + by + c = 0
           c) { #coefficient in line equation ax + by + c = 0
    
    #Function calculates shortest distance between point and line in R^2.
    
    if (!(all(is.finite(c(x, y, a, b, c))))) {
      stop("All inputs to linePtDist must be real numbers.")
    }
    
    return(abs(a * x + b * y + c) / sqrt(a ^ 2 + b ^ 2))
  }


callingWaterfall <-
  function (x, type=c("IC50", "AUC", "AMAX"), intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, name="Drug", plot=FALSE) {
    
    type <- match.arg(type)
    if (any(!is.na(intermediate.fold) & intermediate.fold < 0)) { intermediate.fold <- intermediate.fold[!is.na(intermediate.fold) & intermediate.fold < 0] <- 0 }
    if (is.null(names(x))) { names(x) <- paste("X", 1:length(x), sep=".") }
    
    xx <- x[complete.cases(x)]
    print("ok1")
    switch (type,
            "IC50" = {
              xx <- -log10(xx)
              ylabel <- "-log10(IC50)"
              ## 4 fold difference around IC50 cutoff
              if (length(intermediate.fold) == 3) { intermediate.fold <- intermediate.fold[1] }
              if (intermediate.fold != 0) {
                interfold <- log10(intermediate.fold)
              } else { 
                interfold <- 0 
              }
            },
            "AUC" = {
              ylabel <- "AUC"
              ## 1.2 fold difference around Activity Area cutoff
              if (length(intermediate.fold) == 3) { intermediate.fold <- intermediate.fold[2] }
              interfold <- intermediate.fold
            },
            "AMAX" = {
              ylabel <- "Amax"
              ## 1.2 fold difference around Amax
              if (length(intermediate.fold) == 3) { intermediate.fold <- intermediate.fold[3] }
              interfold <- intermediate.fold
            }
    )
    print("ok2")
    if (length(xx) < 3) {
      tt <- array(NA, dim=length(x), dimnames=list(names(x)))
      if (interfold == 0) {
        tt <- factor(tt, levels=c("resistant", "sensitive"))
      } else {
        tt <- factor(tt, levels=c("resistant", "intermediate", "sensitive"))
      }
      return (tt)
    }
    
    oo <- order(xx, decreasing=TRUE)
    ## test linearity with Pearson correlation
    cc <- stats::cor.test(-xx[oo], 1:length(oo), method="pearson")
    ## line between the two extreme sensitivity values
    dd <- cbind("y"=xx[oo][c(1, length(oo))], "x"=c(1, length(oo)))
    rr <- lm(y ~ x, data=data.frame(dd))
    ## compute distance from sensitivity values and the line between the two extreme sensitivity values
    ddi <- apply(cbind(1:length(oo), xx[oo]), 1, function(x, slope, intercept) {
      return(distancePointLine(x=x[1], y=x[2], a=slope, b= -1, c=intercept))
    }, slope=rr$coefficients[2], intercept=rr$coefficients[1])
    if(cc$estimate > cor.min.linear){
      ## approximately linear waterfall
      returnval <- min(abs(xx[oo] - median(xx[oo])))
      cutoff <- which.min(abs(xx[oo] - median(xx[oo])))
      cutoffn <- names(cutoff)[1]
    } else {
      ## non linear waterfall
      ## identify cutoff as the maximum distance
      returnval <- max(abs(ddi)) 
      cutoff <- which.max(abs(ddi))
      cutoffn <- names(ddi)[cutoff]
    }
    ## identify intermediate sensitivities
    switch (type,
            "IC50" = {
              if (interfold == 0) {
                rang <- c(xx[oo][cutoff], xx[oo][cutoff])
              } else {
                rang <- c(xx[oo][cutoff] - interfold, xx[oo][cutoff] + interfold)
              }
            },
            "AUC" = {
              if (interfold == 0) {
                rang <- c(xx[oo][cutoff], xx[oo][cutoff])
              } else {
                rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
              }
            },
            "AMAX" = {
              if (interfold == 0) {
                rang <- c(xx[oo][cutoff], xx[oo][cutoff])
              } else {
                rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
              }
            }
    )
    
    
    ## check whether range is either min or max
    if (rang[2] >= max(xx)) {
      rang[2] <- sort(unique(xx), decreasing=TRUE)[2]
    }
    if (rang[2] <= min(xx)) {
      rang[2] <- sort(unique(xx), decreasing=FALSE)[2]
    }
    if (rang[1] <= min(xx)) {
      rang[1] <- sort(unique(xx), decreasing=FALSE)[2]
    }
    if (rang[1] >= max(xx)) {
      rang[1] <- sort(unique(xx), decreasing=TRUE)[2]
    }
    
    ## compute calls
    calls <- rep(NA, length(xx))
    names(calls) <- names(xx)
    calls[xx < rang[1]] <- "resistant"
    calls[xx >= rang[2]] <- "sensitive"
    calls[xx >= rang[1] & xx < rang[2]] <- "intermediate"
    
    if (plot) {
      par(mfrow=c(2, 1))
      ccols <- rainbow(4)
      mycol <- rep("grey", length(xx))
      names(mycol) <- names(xx)
      mycol[calls == "sensitive"] <- ccols[2]
      mycol[calls == "intermediate"] <- ccols[3]
      mycol[calls == "resistant"] <- ccols[4]
      mycol[cutoffn] <- ccols[1]
      mypch <- rep(16, length(xx))
      names(mypch) <- names(xx)
      mypch[cutoffn] <- 19
      plot(xx[oo], col=mycol[oo], pch=mypch[oo], ylab=ylabel, main=sprintf("%s\nWaterfall", name))
      points(x=cutoff, y=xx[cutoffn], pch=mypch[cutoffn], col=mycol[cutoffn])
      graphics::abline(a=rr$coefficients[1], b=rr$coefficients[2], lwd=2, col="darkgrey")
      lines(x=c(cutoff, cutoff), y=c(par("usr")[3], xx[cutoffn]), col="red")
      lines(x=c(par("usr")[1], cutoff), y=c(xx[cutoffn], xx[cutoffn]), col="red")
      legend("topright", legend=c(sprintf("resistant (n=%i)", sum(!is.na(calls) & calls == "resistant")), sprintf("intermediate (n=%i)", sum(!is.na(calls) & calls == "intermediate")), sprintf("sensitive (n=%i)", sum(!is.na(calls) & calls == "sensitive")), "cutoff", sprintf("R=%.3g", cc$estimate)), col=c(rev(ccols), NA), pch=c(16, 16, 16, 19, NA), bty="n")
      
      plot(ddi, pch=mypch[oo], col=mycol[oo], ylab="Distance", main=sprintf("%s\n%s", name, "Distance from min--max line"))
      points(x=cutoff, y=ddi[cutoffn], pch=mypch[cutoffn], col=mycol[cutoffn])
      legend("topright", legend=c("resistant", "intermediate", "sensitive", "cutoff"), col=rev(ccols), pch=c(16, 16, 16, 19), bty="n")
    } 
    
    tt <- rep(NA, length(x))
    names(tt) <- names(x)
    tt[names(calls)] <- calls
    if (interfold == 0) {
      tt <- factor(tt, levels=c("resistant", "sensitive"))
    } else {
      tt <- factor(tt, levels=c("resistant", "intermediate", "sensitive"))
    }
    return(returnval)  
  }

computeInteractionMatrix <- function(sensitivityMatrix, measurementType, tradeOff, plot = F){
  interactionMat <- t(apply(sensitivityMatrix, 1, function(x){ callingWaterfall(x, type = measurementType, intermediate.fold = tradeOff, plot = plot)}))
  interactionMat <- apply(interactionMat, 2, function(x){ifelse(is.na(x) | x=="resistant" | x=="intermediate", yes = 0, no = 1)})
  return(interactionMat)
  
}

