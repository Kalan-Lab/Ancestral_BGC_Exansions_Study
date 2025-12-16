library(ggplot2)

dat <- read.table("Order_Codoff_Percentiles.txt", header=T, sep='\t')

pdf("Codoff_Percentile_Histogram.pdf", height=5, width=4)
ggplot(dat, aes(x=Codoff_Discordance_Percentile)) + geom_histogram(color='black', bins=20) + theme_bw() + scale_x_log10()

dev.off()
