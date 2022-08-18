#' Loss method for "plmm" class
#'
#' @param y Observed response vector
#' @param yhat Predicted response vector
#' @export
#' 
#' @examples 
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X), intercept = FALSE)
#' head(loss.plmm(yhat = (fit$SUy), y = admix$y))
#' #TODO (Aug. 18, 2022): ensure that the above choice of 'yhat' is sensible 
 
 
loss.plmm <- function(y, yhat) {
  val <- (y - yhat)^2
  return(val)
}
