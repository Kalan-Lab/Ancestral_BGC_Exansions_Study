library(ggplot2)
library(cowplot)

dat <- read.table("Phylogenetic_Regression_Input.txt", header=T, sep='\t')
# name    gca     bgcome_size     genome_size     cazy_count      protein_count

pdf("Scatterplots_of_Association.pdf", height=3, width=11)

colors <- c('#C00000', '#0B5394', '#BF9000', '#969C90')
names(colors) <- c('Agaricomycetes', 'BGC_Enriched_Pezizomycotina', 'Neocallimastigomycota', 'Other')

g1 <- ggplot(dat, aes(x=genome_size/1000000, y=bgcome_size/1000000, color=clade)) + geom_point(alpha=0.7, show.legend=F) + theme_bw() + xlab("Genome size (Mbp)") + ylab("BGC-ome size (Mbp)") + scale_color_manual(values=colors) + scale_x_log10()
g2 <- ggplot(dat, aes(x=protein_count, y=bgcome_size/1000000, color=clade)) + geom_point(alpha=0.7, show.legend=F) + theme_bw() + xlab("Protein count") + ylab("BGC-ome size (Mbp)")  + scale_color_manual(values=colors)
g3 <- ggplot(dat, aes(x=cazy_count, y=bgcome_size/1000000, color=clade)) + geom_point(alpha=0.7, show.legend=F) + theme_bw() + xlab("CAZy protein count") + ylab("BGC-ome size (Mbp)") + scale_color_manual(values=colors)

plot_grid(g1, g2, g3, nrow=1)
dev.off()
