library(ggplot2)

dat <- read.table("Actinomycetota_SsgB_Plot_Input.txt", header=F, sep='\t')

png("Actinomycetota_Plot.png", height=5, width=2, units='in', res=600)
ggplot(dat, aes(x=reorder(V3, V4), y=V2)) + geom_jitter(color='#9f70b5', alpha=0.1) + geom_boxplot(fill='#CBB8D8', color='#9f70b5', alpha=1.0, outlier.shape = NA ) + theme_classic() + theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1)) + xlab("") + ylab("BGC-ome size")
dev.off()
