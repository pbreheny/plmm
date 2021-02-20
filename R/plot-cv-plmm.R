#' Plot method for cv.plmm class
#'
#' @param x An object of class cv.plmm
#' @param log.l Logical to indicate the plot should be returned on the natural log scale. Defaults to \code{log.l = FALSE}.
#' @param type Type of plot to return. Defaults to "cve."
#' @param selected Logical to indicate which variables should be plotted. Defaults to TRUE.
#' @param vertical.line Logical to indicate whether vertical line should be plotted at the minimum/maximum value. Defaults to TRUE.
#' @param col Color for vertical line, if plotted. Defaults to "red."
#' @param ... Additional arguments.
#' @export

## from cv.ncvreg
plot.cv.plmm <- function(x, log.l=TRUE, type=c("cve", "rsq", "scale", "snr", "all"), selected=TRUE, vertical.line=TRUE, col="red", ...) {
  type <- match.arg(type)
  if (type=="all") {
    plot(x, log.l=log.l, type="cve", selected=selected, ...)
    plot(x, log.l=log.l, type="rsq", selected=selected, ...)
    plot(x, log.l=log.l, type="snr", selected=selected, ...)
    plot(x, log.l=log.l, type="scale", selected=selected, ...)
    return(invisible(NULL))
  }
  l <- x$lambda
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)

  ## Calculate y
  L.cve <- x$cve - x$cvse
  U.cve <- x$cve + x$cvse
  if (type=="cve") {
    y <- x$cve
    L <- L.cve
    U <- U.cve
    ylab <- "Cross-validation error"
  } else if (type=="rsq" | type == "snr") {
    rsq <- pmin(pmax(1 - x$cve/x$null.dev, 0), 1)
    rsql <- pmin(pmax(1 - U.cve/x$null.dev, 0), 1)
    rsqu <- pmin(pmax(1 - L.cve/x$null.dev, 0), 1)
    if (type == "rsq") {
      y <- rsq
      L <- rsql
      U <- rsqu
      ylab <- ~R^2
    } else if(type=="snr") {
      y <- rsq/(1-rsq)
      L <- rsql/(1-rsql)
      U <- rsqu/(1-rsqu)
      ylab <- "Signal-to-noise ratio"
    }
  } else if (type=="scale") {
    y <- sqrt(x$cve)
    L <- sqrt(L.cve)
    U <- sqrt(U.cve)
    ylab <- ~hat(sigma)
  }

  ind <- which(is.finite(l[1:length(x$cve)]))
  ylim <- if (is.null(x$cvse)) range(y[ind]) else range(c(L[ind], U[ind]))
  aind <- intersect(ind, which((U-L)/diff(ylim) > 1e-3))
  plot.args = list(x=l[ind], y=y[ind], ylim=ylim, xlab=xlab, ylab=ylab, type="n", xlim=rev(range(l[ind])), las=1)
  new.args = list(...)
  if (length(new.args)) plot.args[names(new.args)] = new.args
  do.call("plot", plot.args)
  if (vertical.line){
    graphics::abline(v=l[x$min], lty=1, lwd=.5)
    graphics::abline(v=l[x$min1se], lty=2, lwd=.5)
  }
  suppressWarnings(graphics::arrows(x0=l[aind], x1=l[aind], y0=L[aind], y1=U[aind], code=3, angle=90, col="gray80", length=.05))
  graphics::points(l[ind], y[ind], col=col, pch=19, cex=.5)
  if (selected) {
    n.s <- predict.plmm(x$fit, lambda=x$lambda, type="nvars")
    graphics::axis(3, at=l, labels=n.s, tick=FALSE, line=-0.5)
    graphics::mtext("Variables selected", cex=0.8, line=1.5)
  }
}
