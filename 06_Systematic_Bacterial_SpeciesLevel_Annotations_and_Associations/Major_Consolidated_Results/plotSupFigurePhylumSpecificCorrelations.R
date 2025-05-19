library(ggplot2)
library(cowplot)

df <- read.table('../Comprehensive_Info_for_Bacterial_Species_Reps.tsv', header=T, sep='\t')

colors <- c('#9f70b5', '#687ca8', '#669662', '#8fb4b8', '#ba5f7c', '#c9a675')
names(colors) <- c('Actinomycetota', 'Ascomycota', 'Bacillota', 'Cyanobacteriota', 'Myxococcota', 'Pseudomonadota')


png("SupFigure_Correlations.png", height=12, width=12, units='in', res=600)
g1 <- ggplot(df, aes(x=genome_size, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("Genome size") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
g2 <- ggplot(df, aes(x=phage_sum, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("Phage/Mobilome-ome sum") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
g3 <- ggplot(df, aes(x=transposon_count, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("IS element homologs") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
g4 <- ggplot(df, aes(x=cazy_distinct_count, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("Distinct CAZy families") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
g5  <- ggplot(df, aes(x=cazy_total_count, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("Total CAZy homologs") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
g6 <- ggplot(df, aes(x=adherence_count, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("Distinct adherence gene homologs") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
g7 <- ggplot(df, aes(x=oxidative_phospho_count, y=relaxed_bgcome_size, color=phylum)) + theme_bw() + scale_color_manual(values=colors) + facet_wrap(~factor(phylum, levels=c('Bacillota', 'Pseudomonadota', 'Cyanobacteriota', 'Myxococcota', 'Actinomycetota')), scales='free', nrow=1) + xlab("Distinct oxidative phosphorylation enzymes") + ylab("BGC-ome size") + geom_point(show.legend=F, alpha=0.5) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + scale_y_continuous(labels = function(x) format(x, scientific = TRUE))

plot_grid(g1, g2, g3, g4, g5, g6, g7, ncol = 1)

dev.off()
