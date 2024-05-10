library(ggplot2)

dat <- read.table("Gap_Distribution.txt", header=F, sep='\t')

pdf("Gap_Distribution.pdf", height=5, width=5)
ggplot(dat, aes(x=V1)) + geom_histogram() + theme_bw() 
dev.off()
