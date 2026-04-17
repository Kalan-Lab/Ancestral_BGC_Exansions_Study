library(ggplot2)

dat <- read.table("Cyano_Multicellularity_Plot_Input.txt", header=T, sep='\t')

png("Cyano_Plot_V2.png", height=3, width=10, units='in', res=600)
ggplot(dat, aes(x=as.factor(Cyanobacteriota_multicellularity_genes), y=relaxed_bgcome_size, color=Cyanobacteriota_sporulation_genes)) + theme_classic() + theme() + xlab("") + ylab("BGC-ome size") + geom_boxplot(color='#8FB4B8', fill='#8FB4B8', alpha=0.5, outlier.shape=NA) + geom_jitter(show.legend=F, alpha=0.8) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + scale_color_gradient(low='#bcd1d4', high='#557c82')
dev.off()
