library(ggplot2)

dat <- read.table("Overview_Order_Level_with_MedianBGCcounts.with_coloring_info.txt", header=T, sep='\t')

colors <- c('#8FB4B8', '#A070B5', '#BA5E7C')
names(colors) <- c('Cyano', 'Actino', 'Myxo')

pdf("Order_Level_Scatterplot_View.pdf", height=5, width=5)
ggplot(dat, aes(x=(median_bgcome_size_relaxed+0), y=(median_bgc_count_relaxed+0), color=color)) + geom_point(alpha=0.8, size=5, show.legend=F) + theme_bw() + scale_color_manual(values=colors) + geom_hline(yintercept=5, linetype=2, color='red') + geom_vline(xintercept=150000, linetype=2, color='red') + xlab("Median BGC-ome size") + ylab("Median BGC counts") +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

print(dim(dat))
dat.filt <- dat[dat$median_bgcome_size_relaxed <= 150000,]
dat.filt <- dat.filt[dat.filt$median_bgc_count_relaxed <= 5,]
print(dim(dat.filt))
