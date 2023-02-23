#' Print method for summary.cv.plmm objects
#'
#' @param obj An object of class \code{summary.cv.plmm}
#' @param digits The number of digits to use in formatting output 
#' @param ... 
#'
#' @rdname print.summary.cv.plmm 
#' @export
#'
#' @examples
#' cv_fit <- cv.plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' s <- summary.cv.plmm(cv_fit)
#' print.summary.cv.plmm(s)
print.summary.cv.plmm <- function(obj, digits, ...){
  digits <- if (missing(digits)) digits <- c(2, 4, 2, 2, 3) else rep(digits, length.out=5)
  cat(obj$fit$penalty, "-penalized model with n=", nrow(obj$fit$std_X),"\n", sep="")
  cat("At minimum cross-validation error (lambda=", formatC(obj$lambda.min, digits[2], format="f"), "):\n", sep="")
  cat("-------------------------------------------------\n")
  cat("  Nonzero coefficients: ", obj$nvars, "\n", sep="")
  cat("  Cross-validation error (deviance): ", formatC(min(obj$cve), digits[1], format="f"), "\n", sep="")
  cat("  Scale estimate (sigma): ", formatC(sqrt(obj$cve[obj$min]), digits[5], format="f"), "\n", sep="")
}