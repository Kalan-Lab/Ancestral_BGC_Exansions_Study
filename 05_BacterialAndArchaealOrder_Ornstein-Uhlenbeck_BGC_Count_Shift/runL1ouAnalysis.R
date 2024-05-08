library(l1ou)
library(ape)

tree <- read.tree("Bacterial_and_Archaeal_ModelB.rooted.with_innernodes.tre")
dat <- read.table("Overview_Order_Level.txt", header=T, sep='\t')

tree <- chronos(tree)
rownames(dat) <- dat$order

micro <- adjust_data(tree, dat$median_bgcome_size_relaxed)
eModel <- estimate_shift_configuration(micro$tree, micro$Y, nCores=5)
print(eModel)

nEdges <- Nedge(tree) # total number of edges
ew <- rep(1,nEdges) # to set default edge width of 1
ew[eModel$shift.configuration] <- 3 # to widen edges with a shift

pdf('median_bgcome_size_relaxed.pdf', height=10, width=10)
plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
dev.off()

