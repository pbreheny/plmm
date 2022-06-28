#' Calculate scale by the population standard deviation, without centering
#'
#' This function allows you to scale vectors of a matrix by their population standard deviation without centering; we assume our sample is the population.
#' @param X numeric matrix
#' @export
#' 
#' @examples 
#' M <- matrix(rnorm(25), 5, 5)
#' head(M)
#' M_scaled <- scale_varp(M)
#' head(M_scaled)
#' 
#' X_scaled <- scale_varp(admix$X)
#' admix$X[1:5, 1:7]; X_scaled[1:5, 1:7]
scale_varp <- function(x){
  scale(X, center = FALSE, scale = apply(X, 2, function(j) sqrt(mean(j^2, na.rm = TRUE))))
}
