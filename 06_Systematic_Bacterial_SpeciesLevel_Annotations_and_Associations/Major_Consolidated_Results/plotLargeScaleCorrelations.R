library(corrplot)

dat <- read.table("Data_for_Comprehensive_Correlation_Figure.txt", header=T, sep='\t')
dat <- as.matrix(sapply(dat, as.numeric))

pdf("Large_Scale_Correlations.pdf", height=10, width=12)
M=cor(dat,method="s") #create Spearman correlation matrix

corrplot(M, type = 'lower', order = 'hclust', tl.col = 'black',
         cl.ratio = 0.2, tl.srt = 45, col = COL2('PuOr', 10))
dev.off()
