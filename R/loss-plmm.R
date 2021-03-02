#' Loss method for "plmm" class
#'
#' @param y Observed response vector
#' @param yhat Predicted response vector
#' @export
loss.plmm <- function(y, yhat) {
  val <- (y - yhat)^2
  return(val)
}
