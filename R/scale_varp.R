#' Calculate scale by the population standard deviation, without centering
#'
#' This function allows you to scale vectors of a matrix by their population standard deviation without centering; we assume our sample is the population.
#' @param X numeric matrix
#' @keywords internal
#' 
#' @examples 
#' \dontrun{
#' M <- matrix(rnorm(25), 5, 5)
#' head(M)
#' M_scaled <- scale_varp(M)
#' head(M_scaled)
#' 
#' X_scaled <- scale_varp(admix$X)
#' admix$X[1:5, 1:7]; X_scaled[1:5, 1:7]
#' }
scale_varp <- function(X){
  # calculate scaling values
  scale_vals <- apply(X, 2, function(j) sqrt(mean(j^2, na.rm = TRUE)))
  # scale X matrix 
  scaled_X <- scale(X, center = FALSE, scale = scale_vals)
  
  return(list(scale_vals = scale_vals,
              scaled_X = scaled_X))
  
}
