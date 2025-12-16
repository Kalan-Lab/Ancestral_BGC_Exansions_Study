library(ggplot2)

colors <- c('#8FB4B8', '#A070B5', '#BA5E7C', '#858585')
names(colors) <- c('Cyano', 'Actino', 'Myxo', 'Other')

dat <- read.table('biosynrule_account_overlap_psaps_track_modified_for_plotting.txt', header=T, sep='\t')

pdf("Psaps_BioSynRule_Bacterial_Orders_Overlap_Accounted.pdf", height=4, width=4)
ggplot(dat, aes(y=tot_ogs, x=branch_sum, color=group_color)) + geom_point(size=4, alpha=0.8, show.legend=F) + theme_bw() + scale_color_manual(values=colors)
dev.off()
