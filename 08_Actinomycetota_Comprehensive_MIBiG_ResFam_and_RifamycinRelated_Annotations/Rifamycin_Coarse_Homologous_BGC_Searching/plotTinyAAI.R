library(ggplot2)

dat <- read.table('Tiny_AAI_for_Plotting.txt', header=T, sep='\t')

png("Tiny_AAI_for_Plotting.png", height=3, width=10, units='in', res=600)

colors <- c('#CC66E6', '#FF66CC', '#9966FF', '#abaaa9')
names(colors) <- c('Actinomycetia', 'Actinomycetia - Clade-1', 'Actinomycetia - Clade-2', 'Other Actinomycetota')

ggplot(dat, aes(x=AAI, color=ClassGroup, y=Prop_Genes_Found)) + geom_point(size=5, alpha=0.7, show.legend=F) + theme_bw() + scale_color_manual(values=colors) + ylab("") + xlab("") + theme(text = element_text(size=20))

dev.off()
