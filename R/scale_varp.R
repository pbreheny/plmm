#' Calculate scale by the population standard deviation, without centering
#'
#' This function allows you to scale vectors of a matrix by their population standard deviation without centering; we assume our sample is the population.
#' @param x numeric matrix
#' @export
scale_varp <- function(x){
  scale(x, center = FALSE, scale = apply(x, 2, function(y) sqrt(mean(y^2, na.rm = TRUE))))
}
