library(ggplot2)
library(cowplot)

dat <- read.table("Overview_GCA_Level.Filtered.txt", header=T, sep='\t')
# Assembly_ID     Order   Genus   BGC_Count       BGCome_Size     Genome_Size     BGC_Genome_Prop Metallophore_Size       NRPS_Size       PKS_Size        NRPS_or_PKS_Size        Relaxed_BGCome_Size     Relaxed_BGC_Count     Distinct_KOfams Distinct_CAZy   Total_CAZy

g1 <- ggplot(dat, aes(x=Relaxed_BGCome_Size, y=Relaxed_BGC_Count)) + geom_point(alpha=0.4) + theme_bw() + xlab("BGC-ome size") + ylab("BGC count") + geom_smooth(method='lm')
g2 <- ggplot(dat, aes(x=Relaxed_BGCome_Size, y=NRPS_or_PKS_Size)) + geom_point(alpha=0.4)  + theme_bw() + xlab("BGC-ome size") + ylab("Strict NRPS or PKS-ome size") + geom_smooth(method='lm')
g3 <- ggplot(dat, aes(x=Relaxed_BGCome_Size, y=BGC_Genome_Prop)) + geom_point(alpha=0.4) + theme_bw() + xlab("BGC-ome size") + ylab("Proportion of genome as BGCs") + geom_smooth(method='lm')

pdf("Alternate_Metrics_Assessment.pdf", height=3, width=9)
plot_grid(g1, g2, g3, cols=3, rel_widths=c(0.5, 0.5, 0.5))
dev.off()
