library(ggplot2)

dat <- read.table('Plotting_Input.txt', header=T, sep='\t')
# clade_comparison        clade   genome  bgcome_size

pdf("SupFigure_S19.pdf", height=8, width=22)
ggplot(dat, aes(x=clade, y=bgcome_size)) + geom_boxplot() + geom_dotplot(binaxis='y', stackdir='center') + theme_bw() + facet_wrap(~clade_comparison, scales='free') + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + ylab("BGC-ome size")
dev.off()

