#' Calculate the population variance
#'
#' This function allows you to calculate the population variance; we assume our sample is the population.
#' @param x numeric vector
#' @export
#' 
#' @examples 
#' v <- rnorm(5)
#' varp(v)
varp <- function(x) mean((x-mean(x))^2)
