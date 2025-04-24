library(ggplot2)

dat <- read.table("tmp", header=T, sep='\t')
#dat <- read.table("Plotting_Input.txt", header=T, sep='\t')
# order	median_genome_size	median_bgcome_size	median_bgc_count	median_bgcome_size_relaxed	median_bgc_count_relaxed	median_npbgcome_size	median_nrps_count	median_pks_count	median_oxy_count	median_distinct_cazy_count	median_total_cazy_count

colors <- c('#8FB4B8', '#A070B5', '#BA5E7C', '#000000')
names(colors) <- c('Cyano', 'Actino', 'Myxo', 'TMP')

pdf("Order_Level_Scatterplot_View.pdf", height=5, width=5)

ggplot(dat, aes(x=(median_bgcome_size_relaxed+0), y=(median_bgc_count_relaxed+0), color=color)) + geom_point(alpha=0.8, size=5, show.legend=F) + theme_bw() + scale_color_manual(values=colors) + geom_hline(yintercept=5, linetype=2, color='red') + geom_vline(xintercept=150000, linetype=2, color='red') + xlab("Median BGC-ome size") + ylab("Median BGC counts") +  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

print(dim(dat))

dat.filt <- dat[dat$median_bgcome_size <= 150000,]
dat.filt <- dat.filt[dat.filt$median_bgc_count <= 5,]
print(dim(dat.filt))
dev.off()
