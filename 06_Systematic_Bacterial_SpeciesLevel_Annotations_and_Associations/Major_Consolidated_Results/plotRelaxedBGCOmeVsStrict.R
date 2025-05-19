#library(corrplot)
library(ggplot2)

dat <- read.table("Merged_Relaxed_and_Strict_Stats.txt", header=T, sep='\t')
dm <- as.matrix(sapply(dat, as.numeric))

png("Relaxed_vs_Strict.png", height=5, width=5, units='in', res=600)
#M=cor(dm,method="s") #create Spearman correlation matrix

#corrplot(M, type = 'lower', tl.col = 'black',
#         cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))

ggplot(dat, aes(x=BGC_Sum_Relaxed, y=NRPS_or_PKS_Sum_Strict)) + geom_point(alpha=0.4) + geom_abline(slope=1) + theme_bw() + xlab("BGC-ome size") + ylab("Strict NRPS+PKS BGC region size")
dev.off()
