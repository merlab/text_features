# generate supp tables
library(readxl)
library(ggplot2)
library(ggpubr)
library(writexl)
library(tidyr)
library(grid)
library(gridExtra)
library(kableExtra)
library(xtable)
library(pander)
count_decimals <- function(x) {
  # length zero input
  if (length(x) == 0) {
    return(numeric())
  }

  # count decimals
  x_nchr <- x %>%
    abs() %>%
    as.character() %>%
    nchar() %>%
    as.numeric()
  x_int <- floor(x) %>%
    abs() %>%
    nchar()
  x_nchr <- x_nchr - 1 - x_int
  x_nchr[x_nchr < 0] <- 0

  x_nchr
}

formatXLSX <- function(x, valCol = "pearsonCor", pVal = 0.05) {
  inMat <- x[, c("drug", "feature selection method", valCol)] # "feature selection method"
  colnames(inMat) <- c("drug", "Type", "value")
  inMat <- as.data.frame(inMat)
  wideMat <- tidyr::pivot_wider(inMat, names_from = "Type", values_from = "value")
  wideMat <- as.data.frame(wideMat)
  rownames(wideMat) <- wideMat$drug
  wideMat$drug <- NULL
  means <- c()
  sds <- c()
  for (i in colnames(wideMat)) {
    x <- as.vector(unlist(wideMat[, i]))
    means[i] <- mean(x, na.rm = TRUE)
    sds[i] <- sd(x, na.rm = TRUE)
  }
  z <- abs(qnorm(pVal * 8))
  n <- nrow(wideMat)
  sme <- sds / sqrt(n)


  # means <- signif(means, 4)
  # means <- signif(means, 4)
  # dec <- min(count_decimals(means))
  dec <- 3
  ciMat <- data.frame(
    drug = colnames(wideMat),
    mean = (round(means, dec)),
    lowCI = (round(means - (z * sds / sqrt(n)), dec)),
    highCI = (round(means + (z * sds / sqrt(n)), dec)),
    sme = round(sme, dec)
  )
  return(ciMat)
}
h <- 4

methods <- c("pearsonCor", "spearmanCor", "kendallCor", "RMSE", "MSE", "MAE") #

generateTable <- function(inFile, inSheet, oldRes = FALSE, pVal = 0.05) {
  space <- "  "
  tableMat <- data.frame()
  rowOrder <- c(
    "var-100", "var-500", "L1000", "L1000-tm", "cor-500",
    "RFE", "MRMR", "GA",
    "text-mining"
  )
  for (i in methods) {
    raw <- read_xlsx(inFile, sheet = inSheet)
    input <- formatXLSX(raw, valCol = i, pVal = pVal)
    if (oldRes == TRUE && i == "pearsonCor") {
      input <- formatXLSX(raw, valCol = i, pVal = 0.00001)
    }
    ref <- input["text-mining", "mean"]
    input$sig <- (input$highCI <= ref | input$lowCI >= ref)
    input$sig <- ifelse(input$sig == TRUE, "*", "")
    input$text <- paste0(space, input$mean, " (±", input$sme, ")", input$sig, space)
    tableMat[input$drug, i] <- input$tex # t
  }
  ylabs <-
    c(
      "pearsonCor" = "Pearson",
      "spearmanCor" = "Spearman",
      "kendallCor" = "Kendall",
      "RMSE" = "RMSE",
      "MSE" = "MSE",
      "MAE" = "MAE"
    )
  colnames(tableMat) <- ylabs[colnames(tableMat)]
  row_colors <- rep(c("white", "grey90"), length.out = nrow(tableMat))
  coln <- colnames(tableMat)
  newColn <- "Feature selection\nmethod"
  tableMat[, newColn] <- rownames(tableMat)
  tableMat <- tableMat[, c(newColn, coln)]

  tableMat <- tableMat[rev(rowOrder), ]
  tableMat <- na.omit(tableMat)

  grid.table(tableMat,
    rows = NULL,
    theme = ttheme_default(
      core = list(
        bg_params = list(fill = row_colors),
        size = 2
      ),
      rowhead = list(
        # bg_params = list(fontface = "plain") # Set row names font to plain
      ),
      hline = list(col = "grey", lwd = 0.5), # Horizontal lines
      vline = list(col = "grey", lwd = 0.5), # Vertical lines
      # padding = unit(c(4, 2), "mm")
    )
  )
}
pdf("./result/supp_tables.pdf", height = 3.5, width = 12)
f <- "./result/Supplementary-data-2-performance-indexes.xlsx"
# f <- "./mlModelMetricsAllCellLine.xlsx"
generateTable(
  inFile = f, inSheet = "RandomForest perf metrics-test",
  pVal = 0.05,
  oldRes = TRUE
)
grid.newpage()
generateTable(
  inFile = f, inSheet = "ElasticNet perf metrics-test",
  pVal = 0.05
)
grid.newpage()
generateTable(
  inFile = "./result/Supplementary-data-2-performance-indexes.xlsx",
  inSheet = "DeepLearning perf metrics",
  pVal = 0.05
)


dev.off()

print("done")
