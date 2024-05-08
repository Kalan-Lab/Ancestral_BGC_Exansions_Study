library(ggplot2)
library(cowplot)

dat <- read.table("Boxplot_Inputs.txt", header=T, sep='\t')
# GCA     Phylum_or_Clade Complete_BGC_Count      BGCome_Size     Strict_NRPS_or_PKSome_Size      Color   Grid

colors <- c('#F4CCCC', '#c9b967', '#737373')
names(colors) <- c('haploid-dominant', 'polyploid-dominant', 'variable')


dat.pez <- dat[dat$Grid == '1 Pez',]
dat.oth <- dat[dat$Grid == '2 other',]


png("Strict_NRPS_or_PKSome_Size.png", res=600, units='in', height=4, width=12)
g1 <- ggplot(dat.pez, aes(x=Phylum_or_Clade, y=Strict_NRPS_or_PKSome_Size)) + geom_dotplot(show.legend=F, aes(fill=Color), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all') + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + scale_fill_manual(values=colors) + geom_boxplot(color='black', aes(fill=Color), alpha=0.3, outlier.shape=NA, show.legend=F)
g2 <- ggplot(dat.oth, aes(x=Phylum_or_Clade, y=Strict_NRPS_or_PKSome_Size)) + geom_dotplot(aes(fill=Color), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all', show.legend=F) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + scale_fill_manual(values=colors) + geom_boxplot(color='black', aes(fill=Color), alpha=0.3, outlier.shape=NA, show.legend=F)
plot_grid(g1, g2, align='h', ncol=2, rel_widths=c(1,6))
dev.off()


png("Complete_BGC_Counts.png", res=600, units='in', height=4, width=12)
g1 <- ggplot(dat.pez, aes(x=Phylum_or_Clade, y=Complete_BGC_Count)) + geom_dotplot(show.legend=F, aes(fill=Color), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all') + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + scale_fill_manual(values=colors) + geom_boxplot(color='black', aes(fill=Color), alpha=0.3, outlier.shape=NA, show.legend=F)
g2 <- ggplot(dat.oth, aes(x=Phylum_or_Clade, y=Complete_BGC_Count)) + geom_dotplot(aes(fill=Color), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all', show.legend=F) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + scale_fill_manual(values=colors) + geom_boxplot(color='black', aes(fill=Color), alpha=0.3, outlier.shape=NA, show.legend=F)
plot_grid(g1, g2, align='h', ncol=2, rel_widths=c(1,6))
dev.off()


png("BGCome_Sizes.png", res=600, units='in', height=4, width=12)
g1 <- ggplot(dat.pez, aes(x=Phylum_or_Clade, y=BGCome_Size)) + geom_dotplot(show.legend=F, aes(fill=Color), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all') + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + scale_fill_manual(values=colors) + geom_boxplot(color='black', aes(fill=Color), alpha=0.3, outlier.shape=NA, show.legend=F)
g2 <- ggplot(dat.oth, aes(x=Phylum_or_Clade, y=BGCome_Size)) + geom_dotplot(aes(fill=Color), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=0.9, binpositions='all', show.legend=F) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + ylab("") + theme(legend.position='bottom') + scale_fill_manual(values=colors) + geom_boxplot(color='black', aes(fill=Color), alpha=0.3, outlier.shape=NA, show.legend=F)
plot_grid(g1, g2, align='h', ncol=2, rel_widths=c(1,6))
#ggsave(file="BGCome_Sizes.svg", plot=pg, height=4, width=12)
dev.off()
