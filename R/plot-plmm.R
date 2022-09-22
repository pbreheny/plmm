#' Plot method for plmm class
#'
#' @param x An object of class \code{plmm}
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the MCP/SCAD penalty and the ridge, or L2 penalty. \code{alpha=1} is equivalent to MCP/SCAD penalty, while \code{alpha=0} would be equivalent to ridge regression. However, \code{alpha=0} is not supported; alpha may be arbitrarily small, but not exactly 0.
#' @param log.l Logical to indicate the plot should be returned on the natural log scale. Defaults to \code{log.l = FALSE}.
#' @param shade Logical to indicate whether a local nonconvex region should be shaded. Defaults to TRUE.
#' @param col Vector of colors for coefficient lines. 
#' @param ... Additional arguments.
#' @export
#' 
#' @examples 
#' fit <- plmm(admix$X[,1:10], admix$y, nlambda = 10)
#' # for the sake of illustration, I consider only 10 SNPs in the plot 
#' plot(fit)
#' plot(fit, log.l = TRUE)

## from ncvreg
plot.plmm <- function(x, alpha=1, log.l=FALSE, shade=TRUE, col, ...) {
  if (length(x$lambda) == 1) stop("Model was fit with only a single lambda value; there is no path to plot", call.=FALSE)
  YY <- if (length(x$penalty.factor)==nrow(x$beta_vals)) coef.plmm(x) else coef.plmm(x)[-1, , drop=FALSE]
  penalized <- which(x$penalty.factor!=0)
  nonzero <- which(apply(abs(YY), 1, sum)!=0)
  ind <- intersect(penalized, nonzero)
  Y <- YY[ind, , drop=FALSE]
  p <- nrow(Y)
  l <- x$lambda

  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else xlab <- expression(lambda)
  plot.args <- list(x=l, y=1:length(l), ylim=range(Y), xlab=xlab, ylab="", type="n", xlim=rev(range(l)), las=1)
  new.args <- list(...)
  if (length(new.args)) plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  if (!is.element("ylab", names(new.args))) graphics::mtext(expression(hat(beta)), side=2, cex=graphics::par("cex"), line=3, las=1)

  if (shade & !is.null(x$convex.min)) {
    l1 <- l[x$convex.min]
    l2 <- min(l)
    graphics::polygon(x=c(l1,l2,l2,l1), y=c(plot.args$ylim[1], plot.args$ylim[1], plot.args$ylim[2], plot.args$ylim[2]), col="gray85", border=FALSE)
  }

  if (missing(col)) {
    col <- grDevices::hcl(h=seq(15, 375, len=max(4, p+1)), l=60, c=150, alpha=alpha)
    col <- if (p==2) col[c(1,3)] else col[1:p]
  } else {
    col <- col[ind]
  }
  line.args <- list(col=col, lwd=1+2*exp(-p/20), lty=1)
  if (length(new.args)) line.args[names(new.args)] <- new.args
  line.args$x <- l
  line.args$y <- t(Y)
  do.call("matlines", line.args)

  graphics::abline(h=0)
}
