library(ggplot2)

dat <- read.table("Codoff_Plotting_Input.txt", header=T, sep='\t')

colors <- c('#c24a82', '#70103d', '#7a7a7a', '#c27c32')
names(colors) <- c('Diploid+ dominant', 'Haploid dominant', 'Dikaryotic', 'Dikaryotic - HET enriched')

fills <- c('#707070', '#b8b6b6', '#f2f2f2', '#f2f2f2', '#f2f2f2')
names(fills) <- c('Complex multicellularity', 'Some complex multicellularity', 'Yeast', 'Zoosporic', 'Other')

pdf("codoff_discordance_percentile.pdf", height=4.5, width=12)
ggplot(dat, aes(x = reorder(clade, codoff_percentile, FUN = median), y=codoff_percentile, color=ploidy_status, fill=morphology)) + geom_boxplot(show.legend=F) + facet_grid(.~large_clade, space='free', scales='free_x') + theme_bw() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + theme(legend.position='bottom') + xlab("") + ylab("Codoff Disturbance Percentile") + scale_fill_manual(values=fills) + scale_color_manual(values=colors) + scale_y_log10()
dev.off()

