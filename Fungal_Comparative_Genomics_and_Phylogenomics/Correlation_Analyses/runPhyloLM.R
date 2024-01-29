library(phylolm)
library(ape)

tree <- read.tree("Phylo_of_BGC-Enriched_Clade.tre")
dat.full <- read.table("Plotting_Input.txt", header=T, sep='\t')
# gca     bgcome_sum      genome_size     cazy_count      starship_count

tree <- chronos(tree)
rownames(dat.full) <- dat.full$gca

trait_info = dat.full$bgcome_sum/1000000
names(trait_info) <- dat.full$gca

pred_info = dat.full$cazy_count
names(pred_info) <- dat.full$gca

print(pred_info)
print(trait_info)
dat = data.frame(trait=trait_info,pred=pred_info)
fit = phylolm(trait~pred,data=dat,phy=tree,model="lambda")
summary(fit)
