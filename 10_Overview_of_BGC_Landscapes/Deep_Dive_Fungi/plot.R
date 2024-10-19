library(ggradar2)

dat <- read.table("Plotting_Input.txt", header=T, sep='\t')

res.dir <- 'Plots/'

groups <- c('BGC_Enriched_Pezizomycotina', 'Agaricomycetes', 'Basal_Pezizomycotina', 'Mucoromycota', 'Pucciniomycotina', 'Saccharomycotina', 'Taphrinomycotina', 'Tremellomycetes', 'Ustilaginomycotina')
fullscores <- c(75.0, 75.0, 75.0, 75.0, 75.0)
for (i in 1:9) {
	group = groups[i]
	print(group)
	dat.filt <- dat[dat$group == group,]
	dat.filt <- dat.filt[,-1]
	row.names(dat.filt) <- dat.filt$cat
	dat.filt <- dat.filt[,-1]
	out.file <- paste(group, '.png', sep='')
	group = row.names(dat.filt)
	df = cbind(group,dat.filt)
	print(out.file)
	png(file=out.file, units='in', res=600, height=6, width=6)
	print(df)
	#max.val <- max(as.numeric(unlist(as.matrix(dat.filt[,-1]))))
	#print(max.val)
	#max.val <- 1.0
	g1 <- ggradar2(df, gridline.label = seq(0, 75, 25), axis.label.offset=1.18, label.centre.y=F, plot.legend=F, group.line.width = 1, centre.y=0.7, group.colours = c('#8d8f8e', '#333333'), group.fill.colours=c('#8d8f8e', '#333333'),  fullscore=fullscores)# grid.min=0.0, grid.max=100.0)
	print(g1)
	dev.off()
}
