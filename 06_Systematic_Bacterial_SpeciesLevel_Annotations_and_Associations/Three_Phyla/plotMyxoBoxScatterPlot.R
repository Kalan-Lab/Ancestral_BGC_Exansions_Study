library(ggplot2)

dat <- read.table("Myxo_ScatterPlot_Input.with_bins.txt", header=T, sep='\t')

png("Myxococcota_Plot.png", height=3, width=2.5, res=600, units='in')
ggplot(dat, aes(x=reorder(social_count_bin, order), y=relaxed_bgcome_size, color=Myxococcota_sporulation_genes)) + geom_boxplot(color='#ba5f7c', fill='#ba5f7c', alpha=0.5, outlier.shape=NA) + geom_jitter(show.legend=F, alpha=0.8) + xlab("") + ylab("BGC-ome size") + theme_classic() + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + scale_color_gradient(low='#ed95b1', high='#9e415f') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
