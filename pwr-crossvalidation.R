# simple function to randomly partition indices 1:n into k equal sized
# groups (folds)
fold <- function(n, k) {
    samp <- sample(n)
    len <- unname(table(cut(seq.int(n), k)))
    end <- cumsum(len)
    mapply(
        function(start, end, samp) {
            samp[start:end]
        },
        start = c(1, head(end, -1)+1),
        end = cumsum(len),
        MoreArgs = list(samp=sample(n)),
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
}

# cross-validation of pwr
pwr.cv <- function(formula, phy4d, wfun, bwidth, method, holdout) {
    # extract training set from phy4d
    phy.train <- phy4d[setdiff(seq(nTips(phy4d)), holdout)]
    if (missing(bwidth)) {
        bwidth <- get.opt.bw(formula, phy.train, wfun=wfun, method=method)
    }
    # get weights vectors for *all* species, but only keep rows
    # corresponding to training set
    wts.train <- getWts(phy4d, bwidth, wfun)[-holdout,]
    # extract training data from phy4d
    dat.train <- tipData(phy.train)
    # loop over each point and do a weighted least squares regression
    # with weights based on distance
    yhat <- sapply(holdout, function(p) {
        w <- wts.train[,p]
        b <- lm(formula, data=data.frame(dat.train, w), weights=w)
        predict(b, tipData(phy4d)[p,])
    })
    return(data.frame(
        y=tipData(phy4d)[holdout, as.character(formula)[[2]]],
        yhat.pwr=yhat)
        )
}


pwr.cv.slope <- function(formula, phy4d, wfun, bwidth, method, holdout) {
    # extract training set from phy4d
    phy.train <- phy4d[setdiff(seq(nTips(phy4d)), holdout)]
    if (missing(bwidth)) {
        bwidth <- get.opt.bw(formula, phy.train, wfun=wfun, method=method)
    }
    # get weights vectors for *all* species, but only keep rows
    # corresponding to training set
    wts.train <- getWts(phy4d, bwidth, wfun)[-holdout,]
    # extract training data from phy4d
    dat.train <- tipData(phy.train)
    # loop over each point and do a weighted least squares regression
    # with weights based on distance
    yhat <- sapply(holdout, function(p) {
        w <- wts.train[,p]
        b <- lm(formula, data=data.frame(dat.train, w), weights=w)
        b$coefficients["x"]
    })
    return(data.frame(
        slope=tipData(phy4d)[holdout, "true_slope"],
	slope.pwr=yhat)
        )
}




# cross-validation of pgls
pgls.cv <- function(formula, phy4d, holdout) {
    # extract training set from phy4d
    p4.train <- phy4d[setdiff(seq(nTips(phy4d)), holdout)]
    phy.train <- suppressWarnings(as(p4.train, "phylo"))
    # extract training data from phy4d
    dat.train <- tipData(p4.train)
    # loop over each point and do a weighted least squares regression
    # with weights based on distance
    pgls <- do.call(gls, list(model=formula, data=dat.train,
        correlation=corBrownian(phy=phy.train)))
    return(data.frame(
        y=tipData(phy4d)[holdout, as.character(formula)[[2]]],
        yhat.pgls=predict(pgls, tipData(phy4d)[holdout,]))
        )
}


pgls.cv.slope <- function(formula, phy4d, holdout) {
    # extract training set from phy4d
    p4.train <- phy4d[setdiff(seq(nTips(phy4d)), holdout)]
    phy.train <- suppressWarnings(as(p4.train, "phylo"))
    # extract training data from phy4d
    dat.train <- tipData(p4.train)
    # loop over each point and do a weighted least squares regression
    # with weights based on distance
    pgls <- do.call(gls, list(model=formula, data=dat.train,
        correlation=corBrownian(phy=phy.train)))
    return(data.frame(
        slope=tipData(phy4d)[holdout, "true_slope"],
        slope.pgls=pgls$coefficients[[2]])
        )
}





simcv <- function(pglsfit, sim=c("slope", "var1", "var2"), vcv, mc.cores) {

    sim <- match.arg(sim, c("slope", "var1", "var2"))

    # make sure target columns don't already exist
    tipData(arn.p4d)$simSlope <- NULL
    tipData(arn.p4d)$seed.sim <- NULL
    tipData(arn.p4d)$ffd.sim <- NULL

    if (sim=="slope") {
        # extract global intercept (assumed constant in our simulation)
        simInt <- coef(pglsfit)[[1]]

        # use global slope estimate as root value for simulated forward
        # evolution of this "trait", and add to tree data
        simSlope <- rTraitCont(as(arn.p4d, "phylo"), sigma=2,
            root.value=coef(pglsfit)[[2]])
        arn.p4d <- addData(arn.p4d, data.frame(simSlope))

        # generate simulated FFD
        arn.p4d <- addData(arn.p4d, data.frame(ffd.sim=simInt +
            tipData(arn.p4d)$simSlope * tipData(arn.p4d)$seed,
            seed.sim=tipData(arn.p4d)$seed))
    } else if (sim=="var1") {
        require(geiger)
        if (missing(vcv)) {
            vcv <- ic.sigma(arnp, tipData(arn.p4d)[c("FFD", "seed")])
        }
        vars.sim <- sim.char(arnp, vcv)[,,1]
        arn.p4d <- addData(arn.p4d, setNames(data.frame(vars.sim), c("ffd.sim",
            "seed.sim")))
    } else if (sim=="var2") {
        require(phytools)
        tipData(arn.p4d)$seed.sim <- fastBM(arnp, a=coef(gls(seed ~ 1,
            data=tipData(arn.p4d), correlation=corBrownian(phy=arnp))))
        tipData(arn.p4d)$ffd.sim <- c(cbind(1, tipData(arn.p4d)$seed.sim) %*%
            coef(pglsfit)) + fastBM(arnp)
    } 


    # 5-fold cross validation
    k <- fold(nTips(arn.p4d), 5)

    if (missing(mc.cores)) {
        mc.cores <- min(length(k), 16)
    }
    # ...pgls
    arn.pgls.cve <- mclapply(seq_along(k), function(fold) {
        yhat <- pgls.cv(ffd.sim ~ seed.sim, arn.p4d, holdout=k[[fold]])
        sqrt(mean(do.call("-", yhat)^2))
        }, mc.cores=mc.cores)
    # ...pwr
    arn.pwr.cve <- mclapply(seq_along(k), function(fold) {
        yhat <- pwr.cv(ffd.sim ~ seed.sim, arn.p4d, wfun="martins",
            method="L-BFGS-B", holdout=k[[fold]])
        sqrt(mean(do.call("-", yhat)^2))
        }, mc.cores=mc.cores)

    return(c(
        sapply(c(mean.pgls=mean, sd.pgls=sd),
            function(f) f(unlist(arn.pgls.cve))),
        sapply(c(mean.pwr=mean, sd.pwr=sd),
            function(f) f(unlist(arn.pwr.cve)))
    ))

}

#
# procedural code
#
stop("end of function definitions")

set.seed(99)
k <- matrix(sample(nTips(ap4d2)), ncol=5)

library(parallel)

# estimate expected generalization error using 5-fold CV
# ...pwr
pgls.cve <- mclapply(seq(ncol(k)), function(fold) {
    yhat <- pgls.cv(ffd.sc ~ seed.sc, ap4d2, holdout=k[,fold])
    sqrt(mean(do.call("-", yhat)^2))
    })
# ...pgls
pwr.cve <- mclapply(seq(ncol(k)), function(fold) {
    yhat <- pwr.cv(ffd.sc ~ seed.sc, ap4d2, wfun="martins",
        method="L-BFGS-B", holdout=k[, fold])
    sqrt(mean(do.call("-", yhat)^2))
    })

sapply(c(mean=mean, sd=sd), function(f) f(unlist(pgls.cve)))
##     mean       sd 
## 39.82692  2.74358 
sapply(c(mean=mean, sd=sd), function(f) f(unlist(pwr.cve)))
##     mean       sd 
## 29.19123  5.25956 

# estimate expected generalization error using 5-fold CV
# ...pwr
arn.pgls.cve <- mclapply(seq(ncol(k)), function(fold) {
    yhat <- pgls.cv(FFD ~ seed, arn.p4d, holdout=k[fold])
    sqrt(mean(do.call("-", yhat)^2))
    })
# ...pgls
arn.pwr.cve <- mclapply(seq(ncol(k)), function(fold) {
    yhat <- pwr.cv(FFD ~ seed, arn.p4d, wfun="martins",
        method="L-BFGS-B", holdout=k[fold])
    sqrt(mean(do.call("-", yhat)^2))
    })

sapply(c(mean=mean, sd=sd), function(f) f(unlist(arn.pgls.cve)))
##     mean       sd 
## 39.82692  2.74358 
sapply(c(mean=mean, sd=sd), function(f) f(unlist(arn.pwr.cve)))
##     mean       sd 
## 29.19123  5.25956 


# estimate expected generalization error using 10-fold CV

set.seed(99)
k <- fold(nTips(arn.p4d), 10)

# ...pgls
arn.pgls.cvp <- mclapply(seq_along(k), function(fold) {
    yhat <- pgls.cv(FFD ~ seed, arn.p4d, holdout=k[[fold]])
    yhat
    }, mc.cores=min(length(k), 16))
arn.pgls.cve <- lapply(arn.pgls.cvp, function(y) {
    sqrt(mean(do.call("-", y)^2))
    })
sapply(c(mean=mean, sd=sd), function(f) f(unlist(arn.pgls.cve)))
##      mean        sd 
## 20.698476  5.332476 
mean(with(do.call(rbind, arn.pgls.cvp), (y-yhat.pgls)^2))
## [1] 454.0954

# ...pwr
arn.pwr.cvp <- mclapply(seq_along(k), function(fold) {
    yhat <- pwr.cv(FFD ~ seed, arn.p4d, wfun="martins",
        method="L-BFGS-B", holdout=k[[fold]])
    yhat
    }, mc.cores=min(length(k), 16))
arn.pwr.cve <- lapply(arn.pwr.cvp, function(y) {
    sqrt(mean(do.call("-", y)^2))
    })
sapply(c(mean=mean, sd=sd), function(f) f(unlist(arn.pwr.cve)))
##      mean        sd 
## 18.491499  3.691288 
mean(with(do.call(rbind, arn.pwr.cvp), (y-yhat.pwr)^2))
## [1] 353.5059



# estimate expected generalization error using 5-fold CV

set.seed(99)
k <- fold(nTips(arn.p4d), 5)

# ...pgls
arn.pgls.cve <- mclapply(seq_along(k), function(fold) {
    yhat <- pgls.cv(FFD ~ seed, arn.p4d, holdout=k[[fold]])
    sqrt(mean(do.call("-", yhat)^2))
    }, mc.cores=min(length(k), 16))
# ...pwr
arn.pwr.cve <- mclapply(seq_along(k), function(fold) {
    yhat <- pwr.cv(FFD ~ seed, arn.p4d, wfun="martins",
        method="L-BFGS-B", holdout=k[[fold]])
    sqrt(mean(do.call("-", yhat)^2))
    }, mc.cores=min(length(k), 16))

sapply(c(mean=mean, sd=sd), function(f) f(unlist(arn.pgls.cve)))
##      mean        sd 
## 21.097889  2.928884 
# ... in-sample prediction error
sqrt(mean(resid(pgls.arn.brownian)^2))
## [1] 21.32004

sapply(c(mean=mean, sd=sd), function(f) f(unlist(arn.pwr.cve)))
##      mean        sd 
## 18.453470  2.916869 

# ... in-sample prediction error
sqrt(mean((mapply(function(coef, seed) coef["(Intercept)",
    "Estimate"] + coef["seed", "Estimate"] * seed, pwr.arn.martins,
    tipData(arn.p4d)$seed) - tipData(arn.p4d)$FFD)^2))
## [1] 8.879647


cve.sim <- mclapply(1:100, function(i) simcv(pgls.arn.brownian, mc.cores=1),
    mc.cores=25)

system.time(cve.sim <- mclapply(1:100, function(i) simcv(pgls.arn.brownian, mc.cores=1),
    mc.cores=25)
)  
##      user    system   elapsed
## 14191.143    10.849   643.536

system.time(cve.sim <- mclapply(1:100, function(i) simcv(pgls.arn.brownian, mc.cores=1),
    mc.cores=25)
)  


#
# simulate data using fastBM
#

library(phytools)
set.seed(99)
tipData(arn.p4d)$x <- fastBM(arnp, a=coef(gls(seed ~ 1, data=tipData(arn.p4d),
    correlation=corBrownian(phy=arnp))))
tipData(arn.p4d)$y <- c(cbind(1, tipData(arn.p4d)$x) %*% coef(gls(FFD ~ seed,
    data=tipData(arn.p4d), correlation=corBrownian(phy=arnp)))) + fastBM(arnp)

# ...pgls
arnBM.pgls.cve <- mclapply(seq_along(k), function(fold) {
    yhat <- pgls.cv(y ~ x, arn.p4d, holdout=k[[fold]])
    sqrt(mean(do.call("-", yhat)^2))
    }, mc.cores=min(length(k), 16))
# ...pwr
arnBM.pwr.cve <- mclapply(seq_along(k), function(fold) {
    yhat <- pwr.cv(y ~ x, arn.p4d, wfun="martins",
        method="L-BFGS-B", holdout=k[[fold]])
    sqrt(mean(do.call("-", yhat)^2))
    }, mc.cores=min(length(k), 16))

arnBM.pgls.yhat <- do.call(rbind, mclapply(1:length(k), function(fold)
    pgls.cv(y ~ x, arn.p4d, holdout=k[[fold]]), mc.cores=5))
arnBM.pwr.yhat <- do.call(rbind, mclapply(1:length(k), function(fold)
    pwr.cv(y ~ x, arn.p4d, wfun="martins", method="L-BFGS-B",
    holdout=k[[fold]]), mc.cores=5))
