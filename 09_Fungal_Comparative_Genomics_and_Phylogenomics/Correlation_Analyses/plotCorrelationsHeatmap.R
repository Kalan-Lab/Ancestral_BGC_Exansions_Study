library(ggplot2)

dat <- read.table("Spearman_Heatmap_Data.txt", header=T, sep='\t')
# Phylum_Order    Phylum  Factor  Factor_Order    Spearman_Correlation    Pval    Pval_Category
png("Spearman_Heatmap.png", height=6, width=10, units='in', res=600)
ggplot(dat, aes(y=reorder(Phylum, -Phylum_Order), x=reorder(Factor, -Factor_Order), color=Spearman_Correlation)) + geom_point(size=6, aes(alpha=Pval_Category)) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_color_gradient2(low='#f54a45', mid='#c2c0c0', high='#4090f7', limits=c(-1, 1)) + xlab("") + ylab("") 
dev.off()
