library(ggplot2)

# Plot comparison between different feature selection methods
df <- data.frame(
  method = c("Text-Mining", "Auto-HMM-LMF", "STF", "SIRS", "ISIRS", "ENR"),
  cor = c(0.608, 0.52, 0.46, 0.48, 0.48, 0.47)
)

clr <- c("#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#66c2a5")
p <- ggplot(data = df, aes(x = method, y = cor, fill = method), alpha = .15) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = rep(clr[2], 6))
p <- p + theme(
  legend.position = "none", axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
  axis.text = element_text(size = 14), panel.background = element_blank(), axis.line.x = element_line(color = "black"),
  axis.line.y = element_line(color = "black")
) +
  xlab("Feature Selection Method") + ylab("Pearson Correlation") + scale_y_continuous(expand = c(0, 0))
p
pdf("./result/Fig4A_Barplot.pdf", width = 6.5, height = 5)
print(p)
dev.off()
