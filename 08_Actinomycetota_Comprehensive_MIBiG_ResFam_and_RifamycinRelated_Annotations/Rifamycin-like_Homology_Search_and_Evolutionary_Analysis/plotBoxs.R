library(ggplot2)

dat <- read.table('Genus_DOG_Jaccard_Indices.txt', header=T, sep='\t')

pdf('Genus_DOG_Jaccard_Indices.pdf', height=4, width=10)
ggplot(dat, aes(x=reorder(comparison, jaccard_index, FUN = median), y=jaccard_index)) + geom_boxplot() + xlab("") + ylab("") 
dev.off()
