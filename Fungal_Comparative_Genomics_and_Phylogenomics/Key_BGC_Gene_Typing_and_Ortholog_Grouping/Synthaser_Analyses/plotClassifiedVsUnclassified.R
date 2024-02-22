library(ggplot2)

dat <- read.table("Classified_vs_Unclassified_by_Clade_Info.txt", header=T, sep='\t')

colors <- c('#5b5c5b', '#aeb0af')
names(colors) <- c('classified', 'unclassified')

pdf("Classified_vs_Unclassified_by_Clade_Info.pdf", height=5, width=10)
ggplot(dat, aes(x=clade, y=count, fill=category)) + geom_bar(color='black', stat='identity', position='dodge') + theme_bw() + xlab("") + ylab("") + scale_fill_manual(values=colors) + scale_y_log10()
dev.off()
