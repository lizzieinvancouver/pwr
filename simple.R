# Regetz et al R code for generating simple tree/data to illustrate the
# value of PWR over conventional global models like PGLS.
#
# Jim Regetz (NCEAS)

require(ape)
require (grid)

source("pwr-functions.R")
source("pwr-plots.R")

# create simple tree as both phylo (tree) and phylo4d (phy) objects,
# with well-behaved (but slightly jittered) predictor and response
# values corresponding to a deep bifurcation in the slope term
set.seed(99)
tree <- read.tree(text="(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);")
dat <- data.frame(
    x = c(0, 1, 2, 3, 0, 1, 2, 3),
    y = rnorm(8, c(1, 2, 3, 4, 2, 1.5, 1, 0.5), sd=0.25))
phy <- phylo4d(tree, dat)

# run pwr with various distance weighting methods
pwr.b <- getEst(pwr(y ~ x, phy, wfun="brownian"))
pwr.m <- getEst(pwr(y ~ x, phy, wfun="martins"))
pwr.g <- getEst(pwr(y ~ x, phy, wfun="gaussian"))
pwr.G <- getEst(pwr(y ~ x, phy, wfun="Gauss"))

# pwr.m:
##          est         lb          ub
## 1  1.0166685  0.8717198  1.16161721
## 2  1.0236991  0.8546364  1.19276175
## 3  1.0520275  0.8860122  1.21804268
## 4  1.0571060  0.9123327  1.20187938
## 5 -0.3922051 -0.5342139 -0.25019636
## 6 -0.3866605 -0.5499172 -0.22340379
## 7 -0.2510523 -0.4700266 -0.03207803
## 8 -0.2793040 -0.4563995 -0.10220841


# run traditional pgls
pgls <- gls(y ~ x, data=dat, correlation=corBrownian(phy=tree))
# print estimates and CIs for pgls
coef(pgls)
## (Intercept)           x
##   1.4281520   0.3106262
confint(pgls)
##                  2.5 %    97.5 %
## (Intercept) -0.1808439 3.0371479
## x           -0.3287239 0.9499764

# generate scatterplot of data
pdf("simple-data.pdf", width=5, height=5)
par(xpd=TRUE)
plot(tipData(phy)$x, tipData(phy)$y, pch=rep(c(1,16), each=4),
    xlab="Predictor", ylab="Response", bty="l")
text(tipData(phy)$x, tipData(phy)$y, labels=rownames(tipData(phy)),
    pos=3, cex=0.8)
dev.off()

# generate plot of pwr (Martins) and pgls estimates
pdf("simple-pwr.pdf", width=5, height=5)
tp(addData(phy, data.frame(gest=coef(pgls)[2], glb=confint(pgls)[2,1],
    gub=confint(pgls)[2,2], pwr.m))[,c("est", "lb", "ub", "gest", "glb",
    "gub")], show.tip.label=TRUE, lower=-1., upper=1.5, pex=2, aex=3.5)
dev.off()

