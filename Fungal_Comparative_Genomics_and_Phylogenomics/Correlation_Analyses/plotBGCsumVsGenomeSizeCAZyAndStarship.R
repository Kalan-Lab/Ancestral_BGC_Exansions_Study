library(ggplot2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')

png("CAZY.png", height=3, width=3, units='in', res=600)
ggplot(dat, aes(x=cazy_count, y=bgcome_sum)) + geom_point(color='#687CA8') + xlab("Distinct CAZy Homologs") +  ylab("BGC-ome size") + theme_bw() + geom_smooth(color='#687CA8', fill='#687CA8', method='lm')
dev.off()

png("Starship.png", height=3, width=3, units='in', res=600)
ggplot(dat, aes(x=starship_count, y=bgcome_sum)) + geom_point(color='#687CA8') + xlab("Starship Conserved Gene Homologs") +  ylab("BGC-ome size") + theme_bw() + geom_smooth(color='#687CA8', method='lm', fill='#687CA8')
dev.off()

png("GenomeSize.png", height=3.5, width=3, units='in', res=600)
ggplot(dat, aes(x=genome_size, y=bgcome_sum)) + geom_point(color='#687CA8') + xlab("") +  ylab("BGC-ome size") + theme_bw() + geom_smooth(color='#687CA8', method='lm', fill='#687CA8') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


