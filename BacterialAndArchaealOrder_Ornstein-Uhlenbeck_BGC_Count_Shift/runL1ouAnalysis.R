library(l1ou)
library(ape)

tree <- read.tree("figtree_simple.with_innernode_names.tre")
dat <- read.table("Overview_Order_Level.txt", header=T, sep='\t')
#ord   median_genome_size      median_bgc_size median_terpene_size     median_oxy_count        median_cazy_count

tree <- chronos(tree)

mbc <- dat$median_bgc_size
names(mbc) <- dat$order
micro <- adjust_data(tree, mbc)
eModel <- estimate_shift_configuration(micro$tree, micro$Y, nCores=5)
print(eModel)

nEdges <- Nedge(tree) # total number of edges
ew <- rep(1,nEdges) # to set default edge width of 1
ew[eModel$shift.configuration] <- 3 # to widen edges with a shift

pdf('median_bgc_count.pdf', height=10, width=10)
plot(eModel, cex=0.5, label.offset=0.02, edge.width=ew)
dev.off()
