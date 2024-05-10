library(ggplot2)

dat <- read.table('Bargraph_Data.txt', header=T, sep='\t')
# Enzyme  Class   Count

cols = c('#303030', '#303030', '#303030', '#CC66E6')
names(cols) <- c('1. Feature DAHP synthase ferredoxin-like domain (PF18152)', '2. Type I DAHP synthases (TIGR01358)', '3. Type II DAHP synthases (TIGR00034)', '4. RifH homolog from AHBA neighborhood')

pdf("Barplot.pdf", height=5, width=10)
ggplot(dat, aes(x=Class, y=Count*100, fill=Enzyme)) + geom_bar(color='black', stat='identity') + facet_wrap(Enzyme~., scales='free', ncol=1) + scale_fill_manual(values=cols) + coord_flip() + theme_bw() + ylab("Percentage of genomes with enzyme")
dev.off()
