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
gic.plmm <- function(fit, ic=c("bic", "hdbic", "ebic"), SUX, SUy, S, eta){

  #need to deal with situations where dim SUX != dim X because of singularity

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

  c(n, p) %<-% dim(SUX[,-1]) # this assumes there is an intercept column
  eta <- fit$eta
  ll <- logLik_nonnull(fit, SUX, SUy, S, eta)
  j <- predict.plmm(fit, type='nvars')
  jj <- pmin(j, p/2) # dont' give smaller penalties for larger models
  df <- j + 2 # +1 for intercept, +1 for sigma2
  # ts <-  choose(p, j) # no penalty if selecting all bc this value is small
  # limit this to select at most half the vars?

  if (ic == "bic"){
    an <- log(n)
    gam <- 0
  } else if (ic == "hdbic"){
    an <- log(log(n))*log(p)
    gam <- 0
  } else if (ic == 'ebic'){
    an <- log(n)
    gam <- 1
  } else {
    stop("IC not implemented.")
  }

  gic <- (-2)*ll + an*df + 2*gam*(lgamma(p+1) - lgamma(j+1) - lgamma(p-j+1))

  return(list(fit = fit,
              lambda = fit$lambda,
              nzero = df - 2,
              gic = gic,
              lambda.min = fit$lambda[which.min(gic)],
              lambda.min.idx = which.min(gic)))
}
