# Davies et al. R code for generating simple tree/data to evaluate the
# performance of PWR across various tree topologies (stemminess and imbalance)
#
# By Jim Regetz (NCEAS)


require(nlme)
require(phylobase)
require(ape)
require(geiger)
require(apTreeshape)
require(parallel)
require(fields)
require(phytools)

source("pwr-functions.R")
source("pwr-crossvalidation.R")

simdir <- "./sims"

set.seed(99)

NSIMS <- 10000

# simulate trees
trees <- lapply(seq_len(NSIMS), function(i) {
        # generate trees with branch lengths drawn from a random uniform
        # distribution
        tree <- rtree(32, br=runif)
        # make tree ultrametric using NPRS
        tree <- chronopl(tree, lambda=1)
        # transform tree stemminess
        tree <- deltaTree(tree, runif(1, 0.1, 0.3))
        # standardise tree root-to-tip length
        rescaleTree(tree, 1)
    })
class(trees) <- "multiPhylo"

# write trees to file
write.tree(trees, file=file.path(simdir, "simulated-trees.phy"))
trees <- read.tree(file.path(simdir, "simulated-trees.phy"))

treesimcv <- function(tree, mc.cores) {

    # tree imbalance
    Ic <- colless(as.treeshape(tree), norm="yule")
    # tree stemminess
    my.gamma <- gammaStat(tree)

    # make phylo4d
    phy <- phylo4d(tree)

sim<-"slope"
intercept<-0
if (sim=="slope") {
  # evolve slope and derive t2 on the simulated tree to evolve
    tipData(phy)$x <- runif(nTips(phy), 0, 1)
    tipData(phy)$true_slope <- fastBM(tree)
    tipData(phy)$y <- intercept + tipData(phy)$x * tipData(phy)$true_slope
    }

    # 5-fold cross validation
    k <- fold(nTips(phy), 5)

    if (missing(mc.cores)) {
        mc.cores <- min(length(k), 16)
    }
    # ...pgls
    pgls.cve <- mclapply(seq_along(k), function(fold) {
        yhat <- pgls.cv.slope(y ~ x, phy, holdout=k[[fold]])
        sqrt(mean(do.call("-", yhat[c("slope","slope.pgls")])^2))
        }, mc.cores=mc.cores)
    # ...pwr
    pwr.cve <- mclapply(seq_along(k), function(fold) {
        yhat <- pwr.cv.slope(y ~ x, phy, wfun="martins",
            method="L-BFGS-B", holdout=k[[fold]])
        sqrt(mean(do.call("-", yhat[c("slope","slope.pwr")])^2))
        }, mc.cores=mc.cores)

    return(c(
        Ic=Ic,
        gamma=my.gamma,
        sapply(c(mean.pgls=mean, sd.pgls=sd),
            function(f) f(unlist(pgls.cve))),
        sapply(c(mean.pwr=mean, sd.pwr=sd),
            function(f) f(unlist(pwr.cve)))
        ))
}

simtrees.cv <- mclapply(trees, treesimcv, mc.cores=5)

Y<-simtrees.cv
Z <- as.data.frame(do.call("rbind", Y[sapply(Y, class)=="numeric"]))
write.csv(Z, file = "PWR_CV_simulations_slope.csv", row.names = F)


pdf(file="CVE-profiles-10ksims_slopes.pdf", height=8.5, width=11)
par(mfrow=c(1,2))
# pwr CV error
scatter.smooth(Z$Ic, Z$mean.pwr, xlab="Ic", ylab="Err(pwr)",
    col="red", ylim = c(0,2), pch=16, cex=0.5)
scatter.smooth(Z$gamma, Z$mean.pwr, xlab="gamma", ylab="Err(pwr)",
    col="red", ylim = c(0,2), pch=16, cex=0.5)
# difference in CV error (pgls-pwr)
scatter.smooth(Z$Ic, Z$mean.pgls-Z$mean.pwr, xlab="Ic",
    ylab="Err(pgls)-Err(pwr)", col="blue", ylim = c(1.5,-1.5), pch=16, cex=0.5)
abline(h=0, lty=2)
scatter.smooth(Z$gamma, Z$mean.pgls-Z$mean.pwr, xlab="gamma",
    ylab="Err(pgls)-Err(pwr)", col="blue", ylim = c(1.5,-1.5), pch=16, cex=0.5)
abline(h=0, lty=2)
dev.off()

