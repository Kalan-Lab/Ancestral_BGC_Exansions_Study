library(ggplot2)

dat <- read.table("Heatmap_Input.txt", header=T, sep='\t')
# Clade   Clade_Order     Domain_ID       Domain  Value

pdf("Heatmap.pdf", heigh=3, width=10)
ggplot(dat, aes(y=reorder(Clade, -Clade_Order), x=Domain, fill=Value)) + geom_tile() + theme_classic() + scale_fill_gradient(low='#232421', high='#f7b228') + xlab("") + ylab("")  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()
