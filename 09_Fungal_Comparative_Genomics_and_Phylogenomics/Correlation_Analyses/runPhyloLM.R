library(phylolm)
library(ape)

tree <- read.tree("Fungi_Wide_Tree_with_Outgroups_pared.tre")
dat.full <- read.table("Pez_CAZy_to_BGComeSize_Data.txt", header=T, sep='\t')
# GCA     GCA_Name	Clade   Rep     BGCome_Size     Genome_Size     Unique_Proteins CAZy_Total_Count

tree <- chronos(tree)
rownames(dat.full) <- dat.full$GCA_Name

trait_info = dat.full$BGCome_Size/1000000
names(trait_info) <- dat.full$GCA_Name

pred_info = dat.full$Genome_Size/1000000
names(pred_info) <- dat.full$GCA_Name

#print(pred_info)
#print(trait_info)
dat = data.frame(trait=trait_info,pred=pred_info)
fit = phylolm(trait~pred,data=dat,phy=tree,model="lambda")
summary(fit)
