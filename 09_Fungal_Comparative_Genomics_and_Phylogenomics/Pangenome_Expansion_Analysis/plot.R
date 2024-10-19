library(ggplot2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')
dat2 <- read.table("BGCome_Size_Plotting_Input.txt", header=T, sep='\t')

colors <- c('#c24a82', '#70103d', '#7a7a7a', '#c27c32')
names(colors) <- c('Diploid+ dominant', 'Haploid dominant', 'Dikaryotic', 'Dikaryotic - HET enriched')

fills <- c('#707070', '#b8b6b6', '#f2f2f2', '#f2f2f2', '#f2f2f2')
names(fills) <- c('Complex multicellularity', 'Some complex multicellularity', 'Yeast', 'Zoosporic', 'Other')

pdf("genome_fluidity.pdf", height=4.5, width=12)
ggplot(dat, aes(x = reorder(group, genome_fluidity, FUN = median), y=genome_fluidity, color=ploidy_status, fill=morphology)) + geom_boxplot(show.legend=F) + scale_y_log10() + facet_grid(.~large_clade, space='free', scales='free_x') + theme_bw() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + theme(legend.position='bottom') + xlab("") + ylab("Unstandardized genome fluidity") + scale_fill_manual(values=fills) + scale_color_manual(values=colors)
dev.off()

pdf("genome_fluidity_standardized.pdf", height=4.5, width=12)
ggplot(dat, aes(x = reorder(group, genome_fluidity_standardized, FUN = median), y=genome_fluidity_standardized, color=ploidy_status, fill=morphology)) + geom_boxplot(show.legend=F) + scale_y_log10() + facet_grid(.~large_clade, space='free', scales='free_x') + theme_bw() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + theme(legend.position='bottom') + xlab("") + ylab("Phylogenetically standardized\ngenome fluidity") + scale_fill_manual(values=fills) + scale_color_manual(values=colors)
dev.off()

pdf("bgcome_size.pdf", height=4.5, width=12)
ggplot(dat2, aes(x = reorder(clade, bgcome_size, FUN = median), y=bgcome_size, color=ploidy_status, fill=morphology)) + geom_boxplot(show.legend=F) + scale_y_log10() + facet_grid(.~large_clade, space='free', scales='free_x') + theme_bw() + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + theme(legend.position='bottom') + xlab("") + ylab("BGC-ome size (Mbp)") + scale_fill_manual(values=fills) + scale_color_manual(values=colors)
dev.off()

