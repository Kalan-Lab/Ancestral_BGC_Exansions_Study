library(ggplot2)

dat <- read.table("MGE_Overlap_with_Length.txt", header=T, sep='\t')
# gca     bgc     phylum  mge_overlap     length


pdf("MGE_Overlap_with_Length.pdf", height=7, width=3)
ggplot(dat, aes(x=length, fill=mge_overlap)) + geom_histogram(bins=10, color='black',show.legend=F) + theme_bw() + scale_x_log10()  + facet_wrap(~phylum, ncol=1, scales='free') + scale_fill_manual(values=c('#a1a1a1', '#e084a0' ))

dev.off()
