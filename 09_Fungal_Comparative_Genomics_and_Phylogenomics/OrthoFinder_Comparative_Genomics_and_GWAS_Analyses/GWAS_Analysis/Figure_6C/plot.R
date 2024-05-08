library(ggplot2)

dat <- read.table("Scatterplot_Input.txt", header=T, sep='\t')

color <- c('#619917', '#177496', '#e04141', '#a1a1a1', '#8d4ca8', '#993284', '#d6862b', '#37b1de', '#bdb764', '#5c9976', '#a8924f', '#ed0047')
names(color) <- c('Carbohydrate utilization', 'Cytochrome', 'Het', 'Hypothetical / other','Methyltransferase','Peptidase / protease','PKS','Redox','Regulatory','Terpene synthase','Transport', 'OTHER')

pdf("Figure_6C_with_Legend.pdf", height=5, width=4)
# HOG     AF      LRT_Pvalue      Annotation_Category
ggplot(dat, aes(x=AF, y=-log(LRT_Pvalue, 10), color=Annotation_Category)) + geom_point(alpha=0.7) + theme_bw() + scale_color_manual(values=color)
dev.off()

png("Figure_6C.png", height=4, width=4, res=600, units='in')
ggplot(dat, aes(x=AF, y=-log(LRT_Pvalue, 10), color=Annotation_Category)) + geom_point(alpha=0.7, show.legend=F) + theme_bw() + scale_color_manual(values=color) + xlab("HOG Frequency") + ylab("LRT P-value (-log10)")
dev.off()
