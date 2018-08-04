# Regetz et al R code defining PWR plot helper functions.
#
# Jim Regetz (NCEAS)

# simple helper function for rescaling values
rescale <- function(x, lower=NULL, upper=NULL, na.rm=TRUE) {
    if (is.null(lower)) {
        lower <- min(x, na.rm=na.rm)
    }
    if (is.null(upper)) {
        upper <- max(x, na.rm=na.rm)
    }
    (x - lower) / (upper - lower)
}

# modified version of phylobase:::treePlot
tp <- function (phy, show.tip.label=TRUE, show.node.label=FALSE,
    tip.order=NULL, tip.plot.fun="bubbles", plot.at.tip=TRUE,
    edge.color="black", node.color="black", tip.color="black",
    edge.width=1, newpage=TRUE, margins=c(1.1, 1.1, 1.1, 1.1), ...) {
    if (!inherits(phy, "phylo4"))
        stop("treePlot requires a phylo4 or phylo4d object")
    if (!isRooted(phy))
        stop("treePlot function requires a rooted tree.")
    if (newpage)
        grid.newpage()
    type <- "phylogram"
    Nedges <- nEdges(phy)
    Ntips <- nTips(phy)
    if (!is.null(tip.order) && length(tip.order) > 1) {
        if (length(tip.order) != Ntips) {
            stop("tip.order must be the same length as nTips(phy)")
        }
        if (is.numeric(tip.order)) {
            tip.order <- tip.order
        } else {
            if (is.character(tip.order)) {
                tip.order <- as.numeric(names(tipLabels(phy))[
                    match(tip.order, tipLabels(phy))])
            }
        }
        tip.order <- rev(tip.order)
    }
    if (!hasEdgeLength(phy) || type == "cladogram") {
        edgeLength(phy) <- rep(1, Nedges)
    }
    xxyy <- phyloXXYY(phy, tip.order)
    pushViewport(plotViewport(margins=margins))
    pb(type=type, show.node.label=show.node.label,
        rot=0, edge.color=edge.color, node.color=node.color,
        tip.color=tip.color, edge.width=edge.width,
        show.tip.label=show.tip.label, newpage=TRUE, ..., XXYY=xxyy)
    upViewport()
}

# modified version of phylobase:::phylobubbles
pb <- function (type=type, place.tip.label="right",
    show.node.label=show.node.label, show.tip.label=show.tip.label,
    edge.color=edge.color, node.color=node.color, tip.color=tip.color,
    edge.width=edge.width, newpage=TRUE, cex=1, pex=1, aex=1, ..., XXYY,
    square=FALSE, show.estimates=TRUE, lower=NULL, upper=NULL) {

    nVars <- 1
    lab.right <- ifelse(place.tip.label %in% c("right", "both"),
        TRUE, FALSE) && show.tip.label
    lab.left <- ifelse(place.tip.label %in% c("left", "both"),
        TRUE, FALSE) && show.tip.label
    phy <- XXYY$phy
    tmin <- min(tdata(phy, type="tip"), na.rm=TRUE)
    tmax <- max(tdata(phy, type="tip"), na.rm=TRUE)
    pedges <- edges(phy)
    tip.order <- XXYY$torder
    tipdata <- tdata(phy, type="tip")[tip.order, , drop=FALSE]
    dlabwdth <- max(stringWidth(colnames(tipdata))) * 1.2
    if (convertWidth(dlabwdth, "cm", valueOnly=TRUE) < 2) {
        dlabwdth <- unit(2, "cm")
    }
    phyplotlayout <- grid.layout(nrow=2, ncol=2, heights=unit.c(unit(1,
        "null"), dlabwdth), widths=unit(c(1, 1), c("null",
        "null"), list(NULL, NULL)))
    pushViewport(viewport(layout=phyplotlayout, name="phyplotlayout"))
    pushViewport(viewport(layout.pos.row=1:2, layout.pos.col=2,
        height=unit(1, "npc") + convertUnit(dlabwdth, "npc"),
        name="bubbleplots", default.units="native"))
    tys <- XXYY$yy[pedges[, 2] <= nTips(phy)]
    tys <- tys[match(names(tipLabels(phy))[tip.order], XXYY$torder)]
    maxr <- ifelse(ncol(tipdata) > nTips(phy), 1/ncol(tipdata),
        1/nTips(phy))
    tipdataS <- apply(tipdata, 2, function(x) (maxr * x)/max(abs(x),
        na.rm=TRUE))
    if (is.null(lower)) {
        lower <- min(tipdata)
    }
    if (is.null(upper)) {
        upper <- max(tipdata)
    }
    tipdataS2 <- rescale(tipdata, lower, upper)
    if (nVars == 1) {
        xpos <- 0.5
    } else {
        xpos <- seq(0 + maxr + 0.12, 1 - maxr - 0.12, length.out=nVars)
    }
    xrep <- rep(xpos, each=length(tys))
    yrep <- rep(tys, nVars)
    ccol <- ifelse(tipdata < 0, "black", "white")
    naxs <- matrix(xrep, ncol=nVars)
    nays <- matrix(yrep, ncol=nVars)
    dnas <- is.na(tipdataS)
    naxs <- naxs[dnas]
    nays <- nays[dnas]
    tipdataS[is.na(tipdataS)] <- 0 + 0.001
    if (lab.right) {
        tiplabwidth <- max(stringWidth(tipLabels(phy)))
    } else {
        tiplabwidth <- unit(0, "null", NULL)
    }
    bublayout <- grid.layout(nrow=2, ncol=2, widths=unit.c(unit(1,
        "null", NULL), tiplabwidth), heights=unit.c(unit(1,
        "null", NULL), dlabwdth))
    pushViewport(viewport(x=0.5, y=0.5, width=0.95, height=1,
        layout=bublayout, name="bublayout"))
    pushViewport(viewport(name="bubble_plots", layout=bublayout,
        layout.pos.col=1, layout.pos.row=1))

    # plot x-axis
    labs <- pretty(c(lower, upper))
    vals <- rescale(labs, lower, upper)
    labs <- format(labs[0 <= vals & vals <= 1], nsmall=1)
    vals <- vals[0 <= vals & vals <= 1]
    ex <- -0.02 * aex
    grid.segments(x0=min(vals), x1=max(vals), y0=ex, y1=ex)
    grid.segments(x0=vals, x1=vals, y0=ex+0.01, y1=ex,
        gp=gpar(col="black"))
    grid.text(labs, x=vals, y=ex-0.01*aex, gp=gpar(cex=0.5*cex))
    grid.text("Coefficient", x=0.5, y=ex-0.02*aex, gp =
        gpar(cex=0.5*cex))

    # plot interesting results
    grid.segments(x0=tipdataS2$gest, x1=tipdataS2$gest, y0=0,
        y1=1, gp=gpar(col="grey"))
    grid.segments(x0=tipdataS2$glb, x1=tipdataS2$glb, y0=0, y1=1,
        gp=gpar(col="grey", lty="dashed"))
    grid.segments(x0=tipdataS2$gub, x1=tipdataS2$gub, y0=0, y1=1,
        gp=gpar(col="grey", lty="dashed"))
    grid.segments(x0=tipdataS2$lb, x1=tipdataS2$ub, y0=tys, y1=tys,
        gp=gpar(col="red"))
    if (!is.null(tipdataS2$lb.1) & !is.null(tipdataS2$ub.1)) {
        grid.segments(x0=tipdataS2$lb.1, x1=tipdataS2$ub.1,
            y0=tys+0.3*diff(tys[1:2]), y1=tys+0.3*diff(tys[1:2]),
            gp=gpar(col="grey"))
    }
    if (!is.null(tipdataS2$lb.2) & !is.null(tipdataS2$ub.2)) {
        grid.segments(x0=tipdataS2$lb.2, x1=tipdataS2$ub.2,
            y0=tys+0.6*diff(tys[1:2]), y1=tys+0.6*diff(tys[1:2]),
            gp=gpar(col="green"))
    }
    if (!is.null(tipdataS2$simcoef)) {
        grid.points(tipdataS2$simcoef, yrep, pch=1, size =
            unit(0.03*pex, "npc"))

    }
    if (show.estimates) {
        grid.points(tipdataS2$est, yrep, pch=16, size=unit(0.02*pex,
            "npc"), gp=gpar(col="red"))

    }

    if (length(naxs) > 0) {
        grid.points(naxs, nays, pch=4)
    }
    upViewport()
    if (lab.right) {
        pushViewport(viewport(name="bubble_tip_labels", layout=bublayout,
            layout.pos.col=2, layout.pos.row=1))
        tt <- sub("_", " ", tipLabels(phy)[tip.order])
        grid.text(tt, 0.1, tys, just="left", gp=gpar(cex=0.5*cex))
        upViewport()
    }
    pushViewport(viewport(name="bubble_data_labels", layout=bublayout,
        layout.pos.col=1, layout.pos.row=2))
    datalaboffset <- convertUnit(unit(15, "mm"), "npc", valueOnly=TRUE)
    upViewport(3)
    pushViewport(viewport(layout.pos.row=2, layout.pos.col=1,
        name="bubblelegend"))
    yyy <- phylobase:::.bubLegendGrob(tipdata, tipdataS)
    grid.draw(yyy)
    upViewport()
    pushViewport(viewport(layout.pos.row=1, layout.pos.col=1,
        name="tree"))
    plotOneTree(XXYY, "phylogram", show.tip.label=FALSE, show.node.label,
        edge.color, node.color, tip.color, edge.width, rot=0)
    upViewport(2)

}
