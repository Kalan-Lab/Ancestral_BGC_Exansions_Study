library(l1ou)
library(ape)

tree <- read.tree("figtree_simple.with_innernode_names.tre")
dat <- read.table("Median_Counts_for_ASR.txt", header=T, sep='\t')

tree <- chronos(tree)
rownames(dat) <- dat$Order

print(class(tree))
micro <- adjust_data(tree, dat[,3])
eModel <- estimate_shift_configuration(micro$tree, micro$Y, nCores=5)
print(eModel)

nEdges <- Nedge(tree) # total number of edges
ew <- rep(1,nEdges) # to set default edge width of 1
ew[eModel$shift.configuration] <- 3 # to widen edges with a shift

pdf('plot.pdf', height=10, width=10)
plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
dev.off()
