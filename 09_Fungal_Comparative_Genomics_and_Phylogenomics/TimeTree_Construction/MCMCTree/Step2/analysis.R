# Code adapted from Mario Dos Reis' github repo divtime: https://github.com/mariodosreis/divtime
# If you use this or further adapt it, please cite: https://pubmed.ncbi.nlm.nih.gov/31278669/

rm(list=ls()) # clean up the workspace

# ###############################################
# POSTERIOR:
# ###############################################
# read in MCMC trace files
mcmc1 <- read.table("mcmc1.txt", head=TRUE)
mcmc2 <- read.table("mcmc2.txt", head=TRUE)
mcmc3 <- read.table("mcmc3.txt", head=TRUE)

# each data frame contains 15 columns:
# MCMC generation number, 9 node ages (divergence times), 2 mean mutation rates,
# 2 rate drift coefficients, and sample log-likelihood values
names(mcmc1)
# [1] "Gen"      "t_n11"    "t_n12"    "t_n13"    "t_n14"    "t_n15"    "t_n16"    "t_n17"
# [9] "t_n18"    "t_n19"    "mu1"      "mu2"      "sigma2_1" "sigma2_2" "lnL"

# to check for convergence of the MCMC runs, we calculate the posterior
# means of times for each run, and plot them against each other
t.mean1 <- apply(mcmc1[,2:76], 2, mean) * 100
t.mean2 <- apply(mcmc2[,2:76], 2, mean) * 100
t.mean3 <- apply(mcmc3[,2:76], 2, mean) * 100

# good convergence is indicated when the points fall on the y = x line.
pdf('Replicability.pdf', height=20, width=20)
par(mfrow=c(2,2))
# posterior times for run 1 vs run 2:
plot(t.mean1, t.mean2, main="a) Posterior times, r 1 vs. r 2"); abline(0, 1)
plot(t.mean1, t.mean3, main="a) Posterior times, r 1 vs. r 3"); abline(0, 1)
dev.off()

# we can calculate the effective sample sizes (ESS) of the parameters
# (you need to have the coda package installed for this to work)
mean.mcmc <- apply(mcmc1[,-1], 2, mean)
ess.mcmc <- apply(mcmc1[,-1], 2, coda::effectiveSize)
var.mcmc <- apply(mcmc1[,-1], 2, var)
se.mcmc <- sqrt(var.mcmc / ess.mcmc)
dat <- cbind(mean.mcmc, ess.mcmc, var.mcmc, se.mcmc)

write.table(dat, file='run_1_info.txt', quote=F, sep='\t', col.names=F)
