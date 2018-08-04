# Davies et al. R code for Arnell PWR analysis. Run successfully using
# R 2.15.1, with external packages phylobase_0.6.3 and ape_3.0-1.
# ...and R 3.4.4, with external packages subplex_1.5-2, phylobase_0.8.4, nlme_3.1-131.1. ape_5.0 
#
# By Jim Regetz (NCEAS)

require(ape)
require(nlme)
require(phylobase)
require(grid)
source("pwr-functions.R")
source("pwr-plots.R")

datadir <- "./data"

# read data and tree
arnd <- read.delim(file.path(datadir, "arnell.txt"),
    stringsAsFactors=FALSE)
arnp <- read.tree(file.path(datadir, "arnell.phy"))

# merge in data to generate phylo4d object
arn.p4d <- phylo4d(arnp, arnd, label.type="column",
    label.column="species")

# first do global pgls
# ...brownian...
pgls.arn.brownian <- gls(FFD ~ seed, data=tipData(arn.p4d),
    correlation=corBrownian(phy=arnp))
# ...martins...
pgls.arn.martins <- gls(FFD ~ seed, data=tipData(arn.p4d),
    correlation=corMartins(value=1, phy=arnp))

# run pwr
# ...brownian...
pwr.arn.brownian <- pwr(FFD ~ seed, arn.p4d, wfun="brownian")
# ...martins...
bw <- get.opt.bw(FFD ~ seed, arn.p4d, wfun="martins", method="subplex")
## bandwidth = 19.9201684626418
pwr.arn.martins <- pwr(FFD ~ seed, arn.p4d, bw=bw, wfun="martins")

# -- manuscript figure --
# plot pgls and martins pwr for real arnell data
pdf("arnell-pwr-m.pdf", width=8, height=9)
tp(addData(arn.p4d[,0], data.frame(
        gest=coef(pgls.arn.martins)[2],
        glb=confint(pgls.arn.martins)[2,1],
        gub=confint(pgls.arn.martins)[2,2],
        getEst(pwr.arn.martins)
    )),
    show.estimates=TRUE, aex=1.2)
dev.off()
