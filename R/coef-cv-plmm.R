#' Coef method for "cv.plmm" class
#'
#' @param object An object of class "cv.plmm."
#' @param lambda A numeric vector of lambda values.
#' @param which Vector of lambda indices for which coefficients to return. Defaults to lambda index with minimum CVE.
#' @param ... Additional arguments (not used).
#' 
#' @rdname coef.cv.plmm
#' 
#' @export
#' 
#' @examples 
#' \dontrun{
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' head(coef.cv.plmm(cv_fit))
#' }
#' 
#' 


coef.cv.plmm <- function(object, lambda, which = object$min, ...){
  coef.plmm(object$fit_to_return, lambda = lambda, which = which, ...)
}
