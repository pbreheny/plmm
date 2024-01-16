#' General information criterion method of selecting lambda for "plmm" class
#'
#' @param fit An object of class "plmm" where svd_detail has been set to TRUE (the default)
#' @param ic Information criterion that should be used to select lambda. Currently supports BIC or HDBIC. Defaults to BIC.
#' @param rot_X Optional: Rotated design matrix including rotated intercept and unpenalized columns, if present. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @param rot_y Optional: Rotated outcome vector. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @param s Optional: Eigenvalues from similarity matrix used for model fitting. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' std_X <- ncvreg::std(admix$X)
#' wUt <- diag(fit$s)%*%fit$U
#' gic_res <- gic(fit = fit, ic = "bic", rot_X = wUt%*%std_X, rot_y = wUt%*%admix$y, s = fit$s)
#' names(gic_res)
#' range(gic_res$gic, na.rm = TRUE) # NAs will result from monomorphic SNPs


gic <- function(fit, ic=c("bic", "hdbic"), rot_X, rot_y, s){
  UseMethod("gic")
}

#' @export
gic.default <- function(fit, ic=c("bic", "hdbic"), rot_X, rot_y, s){
  stop("This function should be used with an object of class plmm_fit")
}

#' @export
gic.plmm <- function(fit, ic=c("bic", "hdbic", "ebic"), rot_X, rot_y, s){

  # TODO: need to deal with situations where dim rot_X != dim X because of singularity

  n <- p <- NULL

  if (missing(rot_X)){
    if (is.null(fit$rot_X)) stop('Rotated X matrix must be supplied as part of fit object or rot_X argument.')
    rot_X <- fit$rot_X
  }

  if (missing(rot_y)){
    if (is.null(fit$rot_y)) stop('Rotated y vector must be supplied as part of fit object or rot_y argument.')
    rot_y <- fit$rot_y
  }

  if (missing(s)){
    if (is.null(fit$s)) stop('Eigenvalues must be supplied as part of fit object or s argument.')
    s <- fit$s
  }


c(n, p) %<-% dim(rot_X[,-1]) # this assumes there is an intercept column
eta <- fit$eta
ll <- log_lik(eta = eta, rot_y = crossprod(fit$U, fit$y), s = s) # TODO: what difference does it make if this is log_lik_nonnull?
j <- predict.plmm(fit, type='nvars')
jj <- pmin(j, p/2) # dont' give smaller penalties for larger models
df <- j + 2 # +1 for intercept, +1 for sigma2
# ts <-  choose(p, j) # no penalty if selecting all bc this value is small
# TODO: consider -- should I limit this to select at most half the vars?


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
