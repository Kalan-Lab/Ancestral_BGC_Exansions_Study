library(ggplot2)

dat <- read.table("Pez_CAZy_to_BGComeSize_Data.txt", header=T, sep='\t')
# GCA     Clade   Rep     BGCome_Size     CAZy_Total_Count

png("CAZY.png", height=2, width=5, units='in', res=600)
ggplot(dat, aes(x=CAZy_Total_Count, y=BGCome_Size)) + geom_point(color='#687CA8', alpha=0.7) + xlab("") +  ylab("") + theme_bw() + geom_smooth(color='#687CA8', fill='#687CA8', method='lm') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

png("Genome_Size.png", height=2, width=5, units='in', res=600)
ggplot(dat, aes(x=Genome_Size, y=BGCome_Size)) + geom_point(color='#687CA8', alpha=0.7) + xlab("") +  ylab("") + theme_bw() + geom_smooth(color='#687CA8', fill='#687CA8', method='lm') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

png("Unique_Proteins.png", height=2, width=5, units='in', res=600)
ggplot(dat, aes(x=Unique_Proteins, y=BGCome_Size)) + geom_point(color='#687CA8', alpha=0.7) + xlab("") +  ylab("") + theme_bw() + geom_smooth(color='#687CA8', fill='#687CA8', method='lm') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

