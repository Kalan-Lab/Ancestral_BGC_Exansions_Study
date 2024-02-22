library(ggplot2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')
# GCA     Domain  Phylum_or_Clade Complete_BGC_Count      Complete_NRPS_or_PKS_BGC_Count

colors <- c('#08345C', '#157CD9', '#0E5595', '#000000', '#78A3AD')
names(colors) <- c('Ascomycota', 'Basidiomycota', 'Pezizomycotina', 'Zoosporic Fungi', 'Zygomycota')

pdf("BigPic.pdf", height=5, width=10)
ggplot(dat, aes(x=Phylum_or_Clade, y=Complete_NRPS_or_PKS_BGC_Count+1)) + geom_boxplot(color='black') + geom_dotplot(aes(fill=Coloring), stackgroups=TRUE, binaxis='y', stackdir='center', dotsize=1.2) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_y_log10() + xlab("") + ylab("") + scale_fill_manual(values=colors)  + theme(legend.position='bottom') 

dev.off()
