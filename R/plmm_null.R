#' Estimate eta using a null Gaussian penalizedLMM model
#'
#' This function allows you to estimate eta (the signal to noise ratio, or narrow-sense variability) under the assumption of a null model.
#' @param K Estimated or known similarity matrix.
#' @param y Continuous outcome vector. If not supplied, V is treated as known and eta is not estimated.
#' @importFrom zeallot %<-%
#' @export

plmm_null<-function(K, y){
  S <- U <- UU <- NULL
  c(S, U, UU) %<-% svd(K)
  Uy <- crossprod(U, y)
  opt <- stats::optimize(f=plmm_nll, c(0.01, 0.99), Uy=Uy, S=S)
  eta <- opt$minimum
  return(list(S=S, U=U, eta=eta))
}
