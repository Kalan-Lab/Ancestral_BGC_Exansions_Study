library(ggplot2)

dat <- read.table("AnnotCategory_Counts.txt", header=T, sep='\t')

color <- c('#619917', '#177496', '#e04141', '#a1a1a1', '#8d4ca8', '#993284', '#d6862b', '#37b1de', '#bdb764', '#5c9976', '#a8924f')
names(color) <- c('Carbohydrate utilization', 'Cytochrome', 'Het', 'Hypothetical / other','Methyltransferase','Peptidase / protease','PKS','Redox','Regulatory','Terpene synthase','Transport')

pdf("AnnotCat_Piechart.pdf", height=7, width=10)
pie(dat$Count, labels=dat$Annotation_Category, border='black', col=color)
dev.off()
