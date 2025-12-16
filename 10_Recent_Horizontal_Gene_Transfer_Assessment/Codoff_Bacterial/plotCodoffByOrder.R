library(ggplot2)

dat <- read.table("Order_Codoff_Percentiles.txt", header=T, sep='\t')

colors <- c('#8FB4B8', '#A070B5', '#A070B5', '#A070B5', '#BA5E7C', '#9a9c9a')
names(colors) <- c('Cyanobacteriales', 'Streptosporangiales', 'Streptomycetales', 'Mycobacteriales', 'Myxococcales', 'Other')

pdf("Order_Codoff_Percentiles.pdf", height=5, width=15)
ggplot(dat, aes(x=reorder(Order, -Codoff_Discordance_Percentile, FUN=median), y = Codoff_Discordance_Percentile, color=Order_Color, fill=Order_Color)) + geom_boxplot(alpha=0.4) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 1 , hjust=1)) + scale_y_log10() + scale_fill_manual(values=colors) + scale_color_manual(values=colors)

dev.off()
