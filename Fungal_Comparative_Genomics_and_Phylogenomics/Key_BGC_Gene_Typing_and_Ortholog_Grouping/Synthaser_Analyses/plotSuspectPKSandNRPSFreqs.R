library(ggplot2)

dat <- read.table("PKS_and_NRPS_Suspect_Frequencies_by_Clades.txt", header=T, sep='\t')

colors <- c('#A64D79', '#E79540')
names(colors) <- c('NRPS(-like)', 'PKS(-like)')

pdf("PKS_and_NRPS_Suspect_Frequencies_by_Clade.pdf", height=5, width=7)
ggplot(dat, aes(x=clade, y=frequency, fill=type)) + geom_bar(color='black', stat='identity', position='dodge') + theme_bw() + xlab("") + ylab("") + scale_fill_manual(values=colors)
dev.off()
