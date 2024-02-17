#' General information criterion method of selecting lambda for "plmm" class
#'
#' @param fit An object of class "plmm" where svd_detail has been set to TRUE (the default)
#' @param ic Information criterion that should be used to select lambda. Currently supports BIC or HDBIC. Defaults to BIC.
#' @param std_X The standardized design matrix
#' @param U Optional: Left singular eigenvectors of X
#' @param s Optional: Eigenvalues from similarity matrix used for model fitting. If not returned as part of \code{plmm} because \code{returnX == FALSE}, must be supplied explicitly.
#' @param eta Optional: Estimate of variance parameter eta 
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples
#' fit <- plmm(X = admix$X, y = admix$y)
#' gic_res <- gic(fit = fit, ic = "bic", std_X = ncvreg::std(admix$X))
#' names(gic_res)
#' range(gic_res$gic, na.rm = TRUE) # NAs will result from monomorphic SNPs

gic <- function(fit, ic=c("bic", "hdbic", "ebic"), std_X, U, s, eta){
  
  # TODO: need to deal with situations where dim rot_X != dim X because of singularity
  # TODO: think about making this a generic method 
  
  n <- p <- NULL
  
  if (missing(eta)){
    if (is.null(fit$eta)) stop('Eta estimate must be supplied as part of fit object or s argument.')
    eta <- fit$eta
  }
  
  if (missing(s)){
    if (is.null(fit$s)) stop('Eigenvalues must be supplied as part of fit object or s argument.')
    s <- fit$s
  }
  
  if (missing(U)){
    if (is.null(fit$U)) stop('Eigenvectors must be supplied as part of fit object or s argument.')
    U <- fit$U
  }
  
  # rotate X (rot_y is returned by fit)
  w <- (eta * s + (1 - eta))^(-1/2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  rot_X <- wUt %*% cbind(1, std_X)
  # re-standardize rotated rot_X, *without* rescaling the intercept!
  # TODO: make this more efficient; should I add an option for fit() to return stdrot_X?
  stdrot_X_temp <- scale_varp(rot_X[,-1, drop = FALSE])
  stdrot_X_noInt <- stdrot_X_temp$scaled_X
  stdrot_X <- cbind(rot_X[,1, drop = FALSE], stdrot_X_noInt) # re-attach intercept
  
  n <- fit$n
  p <- fit$p
  eta <- fit$eta
  ll <- log_lik(eta = eta, rot_y = fit$rot_y, U = U, s = s, n = n) # TODO: what difference does it make if this is log_lik_nonnull?
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
