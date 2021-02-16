#' Estimate eta using a null Gaussian penalizedLMM model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model.
#' @param X Design matrix.
#' @param y Continuous outcome vector.
#' @importFrom zeallot %<-%
#' @export

plmm_null<-function(X, y){
  d <- U <- UU <- NULL
  V <- tcrossprod(ncvreg::std(X))/ncol(X)
  c(d, U, UU) %<-% svd(V)
  #c(d, U) %<-% methods::as(eigen(K), "list")
  Uy <- crossprod(U, y)
  opt <- stats::optimize(f=plmm_nll, c(0.01, 0.99), Uy=Uy, S=d)
  eta <- opt$minimum
  return(list(d=d, U=U, eta=eta, V = V))
}
