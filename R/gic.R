#' General information criterion method of selecting lambda for "plmm" class
#'
#' @param fit An object of class "plmm."
#' @param ic Information criterion that should be used to select lambda. Currently supports BIC or HDBIC. Defaults to BIC.
#' @param SUX Rotated design matrix including rotated intercept and unpenalized columns, if present. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @param SUy Rotated outcome vector. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @param S Eigenvalues from similarity matrix used for model fitting. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @param eta Estimated $eta$ value from object fit.
#' @importFrom zeallot %<-%
#' @export

gic <- function(fit, ic=c("bic", "hdbic"), SUX, SUy, S, eta){
  UseMethod("gic")
}

gic.default <- function(fit, ic=c("bic", "hdbic"), SUX, SUy, S, eta){
  stop("This function should be used with an object of class plmm_fit")
}

#' @export
gic.plmm <- function(fit, ic=c("bic", "hdbic"), SUX, SUy, S, eta){

  n <- p <- NULL

  if (missing(SUX)){
    if (is.null(fit$SUX)) stop('Rotated X matrix must be supplied as part of fit object or SUX argument.')
    SUX <- fit$SUX
  }

  if (missing(SUy)){
    if (is.null(fit$SUy)) stop('Rotated y vector must be supplied as part of fit object or SUy argument.')
    SUy <- fit$SUy
  }

  if (missing(S)){
    if (is.null(fit$S)) stop('Eigenvalues must be supplied as part of fit object or S argument.')
    S <- fit$S
  }

  eta <- fit$eta
  ll <- plmm_nll_nonnull(fit, SUX, SUy, S, eta)
  df <- predict.plmm(fit, type='nvars') + 2 # +1 for intercept, +1 for sigma2
  c(n, p) %<-% dim(SUX)

  if (ic == "bic"){
    an <- log(n)
  } else if (ic == "hdbic"){
    an <- log(log(n))*log(p)
  } else {
    stop("IC not implemented.")
  }

  gic <- (-2)*ll+an*df

  return(list(fit = fit,
              lambda = fit$lambda,
              nzero = df - 2,
              gic = gic,
              lambda.min = fit$lambda[which.min(gic)],
              lambda.min.idx = which.min(gic)))
}
