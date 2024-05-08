library(ggplot2)

dat <- read.table('All_Select_BGC_PropGenomesWith_Info.txt', header=T, sep='\t')
#ST      ST_Type Clade   Prop_Clade_With Pezi_Prop_Clade_With

pdf('NRPS_and_PKS_Subtyping_Plot.pdf', height=9, width=5)

colors <- c('#A64D79', '#E79540', '#c6715c', '#bdbdbd', '#96CEB2')
names(colors) <- c('NRPS', 'PKS', 'Hybrid', '0. Suspected', 'Terpene')

ggplot(dat, aes(x=reorder(ST, -Order), y=Prop_Clade_With, fill=ST_Type)) + geom_bar(color='black',stat='identity') + facet_wrap(~Clade, ncol=1) + theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position='left') + scale_fill_manual(values=colors) + xlab("") + ylab("")

dev.off()
