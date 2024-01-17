library(ggplot2)

dat <- read.table('Simple_Bar_Input.txt', header=T, sep='\t')
# og      group   cat     count   total_actino_count

colors <- c('#a3a2a2', '#a3a2a2', '#000000', '#000000')
names(colors) <- c('other_other', 'other_actino', 'bgc_other', 'bgc_actino')

pdf("Simple_Orthogroup_Counts_Barplot.pdf", height=5, width=20)
ggplot(dat, aes(x=reorder(og, total_actino_count), y=count, fill=cat)) + geom_bar(stat='identity') + theme_classic() + facet_grid(~group, scales='free_x', space='free') + scale_fill_manual(values=colors) +  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), text=element_text(size=20)) + geom_hline(yintercept=0, linetype=2, color='black')
dev.off()
