library(ggplot2)
dat <- read.table('Domain_Architecture_Counts.txt', header=T, sep='\t')

pdf('Domain_Architecture_Counts.pdf', height=7, width=5)
ggplot(dat, aes(x=reorder(Domain_Architecture, Count), y=Count)) + geom_bar(stat='identity', color='#000000', fill='#3b3a3a') + theme_classic() + facet_wrap(~Group, ncol=1, scales='free')  + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) + scale_y_log10()
dev.off()
