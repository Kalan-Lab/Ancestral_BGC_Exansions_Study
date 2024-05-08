library(ggplot2)

dat <- read.table('Comparative_Analysis_Pseudo.txt', header=T, sep='\t')

colors <- c('#FFA10D', '#CB8CAD', '#3E8E93', '#3881B0', '#B292A3', '#DC728C', '#E41A1C', '#C9992C', '#6A886D')
names(colors) <- c('Nocardioidaceae', 'Nanopelagicaceae', 'Dermabacteraceae', 'Propionibacteriaceae', 'Dermatophilaceae', 'Cellulomonadaceae', 'Actinomycetaceae', 'Micrococcaceae', 'Microbacteriaceae')

png("Comparative_Dstat_Boxplots_Pseudo.png", height=4, width=5, units='in', res=600)
# Family  Focal   Part1   Part2   Stat    Total_Count     AABA_Count      ABAA_Count

ggplot(dat, aes(x=Family, y=Stat, color=Family, fill=Family)) + geom_boxplot(alpha=0.7, show.legend=F) + theme_bw() + scale_fill_manual(values=colors) + geom_hline(yintercept=0.0, color='darkgrey', size=1, linetype=2) + ylim(-0.25,0.5) + theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust=1)) + scale_color_manual(values=colors) + ylab("D.stat")

dev.off()
