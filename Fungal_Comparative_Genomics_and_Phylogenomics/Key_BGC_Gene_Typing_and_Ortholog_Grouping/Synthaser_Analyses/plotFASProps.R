library(ggplot2)

dat <- read.table('Proportion_of_Genomes_with_FAS.txt', header=T, sep='\t')

pdf("Proportion_of_Genomes_with_FAS.pdf", height=5, width=7)
ggplot(dat, aes(x=clade, y=prop_clade_with_fas)) + geom_bar(stat='identity', fill='darkgrey', color='black') + theme_bw() + xlab("") + ylab("")
dev.off()
