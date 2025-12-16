library(ggplot2)
library(cowplot)

dat <- read.table('HICount_GenomeSize_BGComeSize_CAZymeCount_Tallies.txt', header=T, sep='\t')
# gca     species bgcome_size     genome_size     hi_count

png('Aspergillus_HET_vs_BGC-ome_Size.png', height=3.5, width=15, units='in', res=600)
g1 <- ggplot(dat, aes(x=hi_count, y=bgcome_size)) + geom_point(alpha=0.7) + theme_bw() + xlab("Proteins with HET or NACHT") + ylab("") + ggtitle("BGC-ome size")
g2 <- ggplot(dat, aes(x=hi_count, y=bgcome_size/genome_size)) + geom_point(alpha=0.7) + theme_bw() + ylab("") + xlab("Proteins with HET or NACHT") + ggtitle("BGC-ome proportion")
g3 <- ggplot(dat, aes(x=hi_count, y=genome_size)) + geom_point(alpha=0.7) + theme_bw() + xlab("Proteins with HET or NACHT") + ylab("") + ggtitle("Genome size")
g4 <- ggplot(dat, aes(x=hi_count, y=cazy_count)) + geom_point(alpha=0.7) + theme_bw() + xlab("Proteins with HET or NACHT") + ylab("") + ggtitle("CAZy count")
plot_grid(g1, g2, g3, g4, nrow=1)#, g4)
dev.off()
