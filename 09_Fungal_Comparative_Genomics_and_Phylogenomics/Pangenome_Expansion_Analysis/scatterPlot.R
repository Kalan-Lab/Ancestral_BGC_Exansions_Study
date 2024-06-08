library(ggplot2)
library(ggrepel)

dat <- read.table('Group_Metrics.txt', header=T, sep='\t')
# statistic       clade   clade_ploidy    number_of_genomes       avg_ogs total_ogs       total_phylo_distance

colors <- c('#ad4c4c', '#c9b967', '#959695')
names(colors) <- c('Haploid Dominant', 'Diploid+ Dominant', 'Variable')


pdf('Fig7.pdf', height=4, width=7)
ggplot(dat, aes(x=total_phylo_distance, y=total_ogs, color=clade_ploidy)) + geom_point(show.legend=F) + theme_bw() + xlim(0,15) + ylim(0,150000) + scale_x_log10() + scale_y_log10() + scale_color_manual(values=colors) + geom_text_repel(aes(label=clade), segment.curvature = -0.1, segment.angle = 20, show.legend=F) + xlab("Branch length sum for phylogeny built\nfrom largely universal proteins") + ylab("Distinct ortholog groups") 
dev.off()
