library(ggplot2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')

pdf("HetProp_SupFigure.pdf", height=4, width=2)
ggplot(dat, aes(x=group, y=(hi_prop*100.0))) + geom_boxplot() + theme_bw()  + theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) + ylab("Percentage of total proteins\nwith HI domain") + xlab("") 
dev.off()


