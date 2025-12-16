library(ggplot2)

dat <- read.table("LCA_PhylogenyDistance_and_CumulativeGeneLoss.txt", header=T, sep='\t')

#colors <- c('#3881B0', '#3E8E93','#449B75','#4AA858','#56A255','#596A98','#6A886D','#7E6E85','#874F6F','#93549D','#999999','#AC5782','#B05B3C','#B16C29','#B292A3','#B53445','#C66764','#C86456','#C9992C','#CB8CAD','#DC728C','#E1C62F','#E3712B','#E41A1C','#E485B7','#F17EB4','#F9F332','#FF7F00','#FFA10D','#FFC31B','#FFE528')
#names(colors) <- c('Propionibacteriaceae', 'Dermabacteraceae', 'Eggerthellaceae', 'Solirubrobacteraceae', 'Ilumatobacteraceae', 'Jiangellaceae', 'Microbacteriaceae', 'Kineococcaceae','Streptomycetaceae','Demequinaceae','Beutenbergiaceae','Mycobacteriaceae','Bifidobacteriaceae','Pseudonocardiaceae','Dermatophilaceae','Rubrobacteraceae','Streptosporangiaceae','Frankiaceae','Micrococcaceae','Nanopelagicaceae','Cellulomonadaceae','Atopobiaceae','Geodermatophilaceae','Actinomycetaceae','Microtrichaceae','Coriobacteriaceae','Micromonosporaceae','Nakamurellaceae', 'Nocardioidaceae', 'Kribbellaceae', 'Brevibacteriaceae')

colors <- c('#FF66CC', '#9966FF') 
names(colors) <- c('Clade-1', 'Clade-2')

png('LCA_PhyloDist_GeneLoss.png', height=6, width=5, units='in', res=600)
ggplot(dat, aes(x=PhyloDist, y=GeneLoss)) + geom_point(aes(color=Clade), size=3, show.legend=F) + geom_smooth(method='lm', show.legend=F) + theme_bw() + scale_color_manual(values=colors) + theme(legend.position='bottom') + xlab("") + ylab("") + theme(text = element_text(size=20))
dev.off()
