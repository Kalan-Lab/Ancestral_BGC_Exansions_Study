library(ggplot2)

dat <- read.table('Plotting_Input.txt', header=T, sep='\t')

pdf('Plot.pdf', height=7, width=9)

colors <- c('#A64D79', '#E79540', '#bdbdbd', '#96CEB2')
names(colors) <- c('NRPS', 'PKS', 'Multiple/Other', 'Terpene')


ggplot(dat, aes(x=reorder(OG_Name, OG_Order), y=Clade_Proportion_With, fill=BGC_Class)) + geom_bar(stat='identity') + facet_wrap(~Clade, ncol=1) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='left') + scale_fill_manual(values=colors) + xlab("") + ylab("")
#OG      OG_Order        Clade   Clade_Proportion_With   BGC_Class       OG_Name

dev.off()
