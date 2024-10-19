library(ggplot2)
library(cowplot)

dat <- read.table('Multi_Type_Stat_Info.txt', header=T, sep='\t')
# clade   tot_bgc_count   multi_type_proportion
dat2 <- read.table('BGC_Length_Info.txt', header=T, sep='\t')
# clade   bgc_cat bgc_length

colors <- c('#A070B5', '#FF5D5D', '#7F9ED7', '#8FB4B8', '#BA5E7C')
names(colors) <- c('Actinomycetota', 'Agaricomycetes', 'BGC_Enriched_Pezizomycotina', 'Cyanobacteriota', 'Myxococcota')

g1 <- ggplot(dat, aes(x=clade, fill=clade, y=multi_type_proportion)) + theme_classic() + geom_bar(color='black', stat='identity', show.legend=F) + scale_fill_manual(values=colors) + xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
g2 <- ggplot(dat2, aes(x=bgc_cat, fill=clade, y=bgc_length)) + theme_bw() + xlab("") + ylab("") + geom_boxplot(color='black', show.legend=F) + scale_fill_manual(values=colors)

png('SupFigure_BC.png', height=4, width=4, units='in', res=600)
plot_grid(g1)
dev.off()
