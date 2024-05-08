library(ggplot2)
library(cowplot)

dat <- read.table("All_Data_for_Sup_Figure.txt", header=T, sep='\t')
# GCA     Phylum_or_Clade Complete_BGC_Count      BGCome_Size     Strict_NRPS_or_PKSome_Size      Color   Grid

dat.a <- dat[dat$Grid == 'a',]
dat.b <- dat[dat$Grid == 'b',]
dat.c <- dat[dat$Grid == 'c',]

png("Targeted_Analyses_SupFigure.png", res=600, units='in', height=5, width=35)
g1 <- ggplot(dat.a, aes(x=Phylum_or_Clade, y=BGCome_Size)) + geom_dotplot(show.legend=F, stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all') + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + geom_boxplot(color='black',  alpha=0.3, outlier.shape=NA, show.legend=F)

g2 <- ggplot(dat.b, aes(x=Phylum_or_Clade, y=BGCome_Size)) + geom_dotplot(show.legend=F, stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all') + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + geom_boxplot(color='black', alpha=0.3, outlier.shape=NA, show.legend=F)


g3 <- ggplot(dat.c, aes(x=Phylum_or_Clade, y=BGCome_Size)) + geom_dotplot(show.legend=F, stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all') + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + geom_boxplot(color='black', alpha=0.3, outlier.shape=NA, show.legend=F)


plot_grid(g1, g2, g3, align='h', ncol=3, rel_widths=c(3,3,3))
dev.off()

