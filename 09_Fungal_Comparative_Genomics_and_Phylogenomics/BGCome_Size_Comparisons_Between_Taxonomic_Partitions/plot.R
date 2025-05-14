library(ggplot2)

dat <- read.table('Plotting_Input.txt', header=T, sep='\t')
# clade_comparison        clade   genome  bgcome_size

pdf("ExtDataFig_7.pdf", height=5, width=20)
ggplot(dat, aes(x=clade, y=bgcome_size)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5, alpha=0.7) + theme_bw() + facet_grid(.~clade_comparison, space='free', scales='free') + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + ylab("BGC-ome size")
dev.off()

