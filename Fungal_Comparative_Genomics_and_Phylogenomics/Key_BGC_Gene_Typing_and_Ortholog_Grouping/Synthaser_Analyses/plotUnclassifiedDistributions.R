library(ggplot2)

dat <- read.table("Unclassified_Breakdown.txt", header=T, sep='\t')

pdf("Unclassified_Breakdown.pdf", height=5, width=10)
ggplot(dat, aes(x=clade, y=count, fill=category)) + geom_bar(stat='identity', color='black', position='fill') + theme_bw() + xlab("") + ylab("") + scale_fill_brewer(palette='Spectral')
dev.off()
