library(ggplot2)

dat <- read.table("Overview_of_BGCOG_Conservation_and_Type.txt", header=F, sep='\t')

png("Barplot.png", height=6, width=7, units='in', res=600)
ggplot(dat, aes(x=V1, y=V3, fill=V2)) + geom_bar(color='black', stat='identity') + theme_classic() + scale_fill_brewer(palette='Reds') + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
