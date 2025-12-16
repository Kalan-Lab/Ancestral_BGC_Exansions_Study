library(ggplot2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')
# domain_evalue   domain_score    meets_ga

pdf("Score_and_Evalue_Histograms.pdf")

ggplot(dat, aes(x=domain_evalue, fill=meets_ga)) + geom_histogram() + theme_bw() + scale_x_log10()
ggplot(dat, aes(x=domain_score, fill=meets_ga)) + geom_histogram() + theme_bw()
dev.off()
