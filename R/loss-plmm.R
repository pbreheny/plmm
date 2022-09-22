#' Loss method for "plmm" class
#'
#' @param y Observed response vector
#' @param yhat Predicted response vector
#' @export
#' 
#' @examples 
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' yhat <- predict.plmm(object = fit, newX = admix$X, type = 'response', lambda = 0.05)
#' head(loss.plmm(yhat = yhat, y = admix$y))

 
 
loss.plmm <- function(y, yhat) {
  val <- (y - yhat)^2
  return(val)
}
