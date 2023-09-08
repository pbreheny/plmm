#' Calculate a relatedness matrix
#'
#' This function allows you to generate an n by n genetic relatedness matrix. If a numeric matrix is supplied, the RRM (Hayes, 2009) is used
#' and is computed XX'/p, where X is standardized. 
#' @param X A numeric matrix of genotypes (from *fully-imputed* data)
#' @export
#' 
#' @examples 
#' RRM <- relatedness_mat(X = admix$X)
relatedness_mat <- function(X){
  # NB: X is standardized as part of the RRM calculation
  rrm <- tcrossprod(ncvreg::std(X))/ncol(X)
  return(rrm)
}

