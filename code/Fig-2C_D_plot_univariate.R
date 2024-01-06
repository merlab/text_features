set.seed(7)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(rstatix)

get_univar_cor_values <- function() {
  dirlc <- "data/Univariate-COR/CCLE"
  fl <- list.files(dirlc, full.names = T)

  fdf <- data.frame(t(sapply(fl, function(i) {
    flx <- rev(strsplit(i, "/")[[1]])[1]
    v <- strsplit(flx, "_")[[1]]
    v
  })))
  colnames(fdf) <- c("drug", "type")
  fdf$file <- rownames(fdf)
  rownames(fdf) <- NULL
  fdf$type <- gsub("\\.rds", "", fdf$type)

  dfl <- list()
  for (dr in unique(fdf$drug))
  {
    tm <- readRDS(fdf$file[fdf$drug == dr & fdf$type == "tm"])
    ntm <- readRDS(fdf$file[fdf$drug == dr & fdf$type == "ntm"])
    dfl[[dr]] <- data.frame(
      Correlation = c(tm, ntm), Type = c(
        rep("Text", length(tm)),
        rep("Other", length(ntm))
      ),
      drug = dr
    )
  }
  return(dfl)
}

createDF <- function(dfl, drugOrd = NA) {
  df <- do.call(rbind.data.frame, dfl)
  df$Type <- factor(df$Type, levels = c("Text", "Other"))

  if (is.na(drugOrd)) {
    mdiff <- sapply(dfl, function(i) {
      mean(i$Correlation[i$Type == "Text"]) - mean(i$Correlation[i$Type == "Other"])
    })
    mdiff <- sort(mdiff, decreasing = F)
    drugOrd <- names(mdiff)
    df$drug <- factor(df$drug, levels = drugOrd)
  }

  mdf <- data.frame(t(sapply(levels(df$drug), function(d) {
    c(
      "Text" = mean(df[df$drug == d & df$Type == "Text", "Correlation"]),
      "Other" = mean(df[df$drug == d & df$Type == "Other", "Correlation"])
    )
  })))

  df$Type <- factor(ifelse(df$Type == "Text", "Text mining", "Other"),
    levels = c("Text mining", "Other")
  )
  colnames(mdf)[1] <- "Text mining"

  return(list(df = df, mean = mdf))
}

## ==============================
## ==============================

dfl <- get_univar_cor_values()
dfallx <- createDF(dfl)

dfall <- dfallx$df
mdx <- dfallx$mean
mdx$drug <- rownames(mdx)
md <- reshape2::melt(mdx, id.vars = "drug", variable.name = "Type")

## =========
md$Type <- factor(as.character(md$Type), levels = c("Text mining", "Other"))

lnc <- "#bababa"
plt <- ggboxplot(md, x = "Type", y = "value", fill = "Type", width = 0.25) +
  geom_line(aes(group = drug), color = lnc, alpha = 0.35) +
  geom_point(aes(fill = Type, group = drug), size = 1.5, shape = 21, color = "white")

# plt + stat_compare_means()

tsr <- wilcox_test(md, value ~ Type)
tsig <- add_significance(tsr)
stat.test <- add_xy_position(tsig, x = "Type")
stat.test <- data.frame(stat.test)
stat.test$ptxt <- paste0("p=", stat.test$p)
plt <- plt + stat_pvalue_manual(stat.test,
  label = "ptxt", size = 2.75,
  tip.length = 0.01, fontface = "italic"
) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.04)))

txtCol <- setNames(
  c( # "#404040", "#ca0020") #"#374E55FF", "#B24745FF"
    "#CD534CFF", "#868686FF"
  ),
  c("Text mining", "Other")
)

plt <- plt + scale_fill_manual(values = txtCol)
plt <- plt + ylab("Mean correlation")

## -------------------------------------------

df <- dfall[dfall$drug %in% rev(levels(dfall$drug))[1:5], ]
rownames(df) <- paste0("r", 1:nrow(df))
df$drug <- droplevels(df$drug)
dodge <- position_dodge(width = 0.8)
## -------------------------

bgCol <- c("#D1E5EB", "#ffffff", "#FFF6E2", "#FAF2E4")[c(1, 2)]
fillCol <- setNames(
  rep_len(bgCol, length.out = length(levels(df$drug))),
  levels(df$drug)
)
drn <- 1:length(fillCol)

pltDr <- ggplot(df) +
  geom_rect(
    data = df[drn, ], xmin = drn - 0.5, xmax = drn + 0.5,
    ymin = -Inf, ymax = Inf, fill = fillCol, alpha = 0.2
  )
pltDr <- pltDr + geom_violin(aes(x = drug, y = Correlation, fill = Type),
  position = dodge,
  alpha = 0.65, size = 0.5, color = "#f0f0f0", trim = T
)
pltDr <- pltDr + geom_boxplot(aes(x = drug, y = Correlation, fill = Type),
  position = dodge,
  notch = F, outlier.size = -1, color = "#252525",
  alpha = 1, lwd = 0.18, width = 0.1
)
pltDr <- pltDr + scale_fill_manual(values = txtCol)
pltDr <- pltDr + scale_x_discrete(expand = c(0, 0)) # , limits = c(0, NA))
# pltDr <- pltDr + scale_y_continuous(expand = c(0, 0), limits = c(0, NA))
pltDr <- pltDr + xlab(NULL)
#### ------------------
stat.test <- df %>%
  group_by(drug) %>%
  wilcox_test(Correlation ~ Type) %>%
  # adjust_pvalue(method = "bonferroni") %>%
  add_significance("p")
stat.test <- stat.test %>% add_xy_position(x = "drug", dodge = 0.8)
stat.test <- data.frame(stat.test)
stat.test$ptxt <- paste0("p=", stat.test$p)
## ---to make all p-value in one line
stat.test$y.position <- max(stat.test$y.position)

pltDr <- pltDr + stat_pvalue_manual(stat.test,
  label = "ptxt", size = 2.0,
  tip.length = 0.01, fontface = "italic", coord.flip = TRUE
) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.001, 0.04)))

pltDr <- pltDr + coord_flip()


#### ------------------
plt <- plt + xlab("") + theme_classic2()
pltDr <- pltDr + theme_classic2()

png("result/Fig-1C_D_univar.png", width = 7.0, height = 6, units = "in", res = 500)
print(ggarrange(plt, pltDr,
  ncol = 2, nrow = 1,
  common.legend = T, legend = c("top", "bottom", "left", "right", "none")[2],
  widths = c(0.3, 0.7)
))
dev.off()
