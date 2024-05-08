library(ggplot2)

dat <- read.table("Spearman_Heatmap_Data_RelaxedBGCOmeSize.txt", header=T, sep='\t')
# Phylum_Order    Phylum  Factor  Factor_Order    Spearman_Correlation    Pval    Pval_Category
png("Spearman_Heatmap.png", height=10, width=8, units='in', res=600)
ggplot(dat, aes(x=reorder(Phylum, -Phylum_Order), y=reorder(Factor, Factor_Order), color=Spearman_Correlation)) + geom_point(size=7, aes(alpha=Pval_Category), show.legend=T) + theme_minimal() + scale_x_discrete(position='top') + theme(axis.text.x = element_text(angle = 45, hjust=0, vjust=0)) + scale_color_gradient2(low='#f54a45', mid='#c2c0c0', high='#4090f7', limits=c(-1, 1)) + xlab("") + ylab("") 
dev.off()

dat <- read.table("Spearman_Heatmap_Data_StrictNRPSandPKSOmeSize.txt", header=T, sep='\t')
# Phylum_Order    Phylum  Factor  Factor_Order    Spearman_Correlation    Pval    Pval_Category
png("Spearman_Heatmap_Strict.png", height=6, width=12, units='in', res=600)
ggplot(dat, aes(y=reorder(Phylum, -Phylum_Order), x=reorder(Factor, -Factor_Order), color=Spearman_Correlation)) + geom_point(size=8, aes(alpha=Pval_Category)) + theme_minimal() + scale_x_discrete(position='top') + theme(axis.text.x = element_text(angle = 45, hjust=0, vjust=0)) + scale_color_gradient2(low='#f54a45', mid='#c2c0c0', high='#4090f7', limits=c(-1, 1)) + xlab("") + ylab("")
dev.off()

