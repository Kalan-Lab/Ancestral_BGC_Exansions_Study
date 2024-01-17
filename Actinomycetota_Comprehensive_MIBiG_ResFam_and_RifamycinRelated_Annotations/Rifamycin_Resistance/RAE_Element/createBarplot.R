library(ggplot2)

dat <- read.table('RAE_Taxonomic_Hit_Distribution.txt', header=T, sep='\t')
# Family  Class   InPhylo ProportionWith  TotalGenomes

png("RAE_Taxonomic_Hit_Distribution.png", height=5, width=15, units='in', res=600)
dat <- dat[dat$TotalGenomes >= 20,]
colors <- c('#CC66E6', '#FF66CC', '#9966FF', '#abaaa9', '#000000', '#FFFFFF')
names(colors) <- c('2Actinomycetia', '3Actinomycetia - Clade-1', '4Actinomycetia - Clade-2', '5Other Actinomycetota', '7With', '6Without')

ggplot(dat, aes(x=reorder(Family, ProportionWith), y=Proportion, fill=WithOrWithout)) + geom_bar(stat='identity', color='white') + theme_classic()  + scale_fill_manual(values=colors)+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("")
dev.off()
