# Regetz et al R code defining core PWR functions.
#
# Jim Regetz (NCEAS)

require(nlme)
require(phylobase)
require(subplex)

# top-level function to run pwr
pwr <- function(formula, phy4d, which.tips, bw, wfun, verbose=FALSE) {
    if (missing(which.tips)) {
        which.tips <- seq.int(nTips(phy4d))
    }
    if (missing(bw)) {
        if (verbose) cat("Computing optimal bandwidth...")
        bw <- get.opt.bw(formula, phy4d, wfun)
        if (verbose) cat("done.\n")
    }
    if (verbose) cat("Using bandwidth of", bw, "\n",
        "Generating weights matrix\n")
    wts <- getWts(phy4d, bw, wfun)
    if (verbose) cat("Running PWR...\n")
    if (verbose) pb <- txtProgressBar(min=0, max=nTips(phy4d), style=3)
    ans <- lapply(which.tips, function(i) {
        if (verbose) setTxtProgressBar(pb, i)
        pwr.wts(formula, phy4d, wts[,i])
    })
    if (verbose) close(pb)
    attr(ans, "weights") <- wts
    return(ans)
}

# lower-level function to run PWR given weight function and bandwidth
pwr.wfn <- function(formula, phy4d, wfun, bwidth, holdout=FALSE) {
    wts <- getWts(phy4d, bwidth, wfun)
    bs <- list()
    res <- c()
    yhat <- c()
    # loop over each point and do a weighted least squares regression
    # with weights based on distance
    for (p in 1:nTips(phy4d)) {
        w <- wts[,p]
        if (holdout) {
            w[p] <- 0
        }
        b <- lm(formula, data=data.frame(tipData(phy4d), w), weights=w)
        bs[[p]] <- coef(b)
        res[p] <- resid(b)[p]
        yhat[p] <- fitted(b)[p]
    }
    coef <- do.call("rbind", bs)
    return(list(coef=coef, resid=res, fitted=yhat))
}


# lower-level function to run PWR given explicit weight values
pwr.wts <- function(formula, phy4d, wts, verbose=FALSE, plot=FALSE) {
    tree <- suppressWarnings(as(phy4d, "phylo"))
    dat <- data.frame(tipData(phy4d), wts=wts)
    # do straight pwr
    pwr.fit <- lm(formula, data=dat, weights=wts)
    if (verbose) {
        cat("\nPWR confidence intervals:\n")
        print(confint(pwr.fit))
    }
    return(coef(summary(pwr.fit)))
}

# optimal bandwidth finder
get.opt.bw <- function(formula, phy4d, wfun, interval, method="subplex",
    verbose=FALSE) {

    if (wfun == "brownian") {
        if (verbose) {
            message("no bandwidth for brownian distance; returning NULL")
        }
        return(NULL)
    }
    # set bounds on possible bandwidth values - these work fine for MS
    # purposes but may need to be revisited for generality
    if (missing(interval)) {
        dist <- cophenetic(suppressWarnings(as(phy4d, "phylo")))
        if (wfun %in% c("gaussian", "Gauss")) {
            lo <- 0.01
            hi <- max(dist)/1.5
        } else if (wfun == "exponential") {
            lo <- 0.01
            hi <- 20
        } else if (wfun == "martins") {
            lo <- 0
            hi <- -log(1e-6)/min(dist[lower.tri(dist)])
        }
        interval <- c(lo, hi)
    }
    if (verbose) message("-- Running ", method, " --")
    if (method=="subplex") {
        optfn <- function(logbwidth) {
            sum(pwr.wfn(formula, phy4d, wfun, exp(logbwidth), TRUE)$resid^2)
        }
        runtime <- system.time(res <- subplex(-1, optfn))
        if (res$convergence!=0) {
           warning(paste("bandwidth optimization problem (code ",
               res$convergence, ")", sep=""))
        }
        opt.bw <- exp(res$par)
    } else if (method=="optimize") {
        optfn <- function(bwidth) {
            sum(pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)$resid^2)
        }
        runtime <- system.time(res <- optimize(optfn, interval=interval))
        opt.bw <- res$minimum
    } else if (method=="L-BFGS-B") {
        optfn <- function(bwidth) {
            sum(pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)$resid^2)
            #vals <- pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)
            #if (any(is.na(vals$coef[, "x"]))) Inf else sum(vals$resid^2)
        }
        runtime <- system.time(res <- optim(1, optfn, lower=0,
            method="L-BFGS-B"))
        opt.bw <- res$par
    } else if (method=="nlm") {
        optfn <- function(bwidth) {
            sum(pwr.wfn(formula, phy4d, wfun, bwidth, TRUE)$resid^2)
        }
        runtime <- system.time(res <- nlm(optfn, 0))
        opt.bw <- res$minimum
    } else {
        stop("invalid optimization algorithm")
    }
    if (verbose) {
        print(runtime)
        message("bandwidth = ", opt.bw)
    }
    opt.bw
}

# create full matrix of weights
getWts <- function(phy4d, bw, wfun) {
    data <- tipData(phy4d)
    tree <- suppressWarnings(as(phy4d, "phylo"))
    dist <- cophenetic(tree)
    switch(wfun,
        gaussian = {
            sqd <- sqrt(dist)
            dnorm(t(t(sqd) / apply(sqd, 2, sd)) / bw)
        },
        exponential = {
            exp(-dist/bw)
        },
        Gauss = {
            exp((-0.5) * ((dist^2)/(bw^2)))
        },
        brownian = {
            corMatrix(Initialize(corBrownian(phy=tree), data))
        },
        martins = {
            corMatrix(Initialize(corMartins(bw, phy=tree), data))
        },
        stop("invalid distance weighting scheme")
    )
}

# extract confidence intervals (handling missingness as needed)
getEst <- function(pwrres) {
    est <- sapply(seq_along(pwrres), function(i) {
        if (nrow(pwrres[[i]])==2) {
            pwrres[[i]][2,1]
        } else NA
    })
    lb <- sapply(seq_along(pwrres), function(i) {
        if (nrow(pwrres[[i]])==2) {
            pwrres[[i]][2,1] - 1.96*pwrres[[i]][2,2]
        } else NA
    })
    ub <- sapply(seq_along(pwrres), function(i) {
        if (nrow(pwrres[[i]])==2) {
            pwrres[[i]][2,1] + 1.96*pwrres[[i]][2,2]
        } else NA
    })
    data.frame(est, lb, ub)
}
