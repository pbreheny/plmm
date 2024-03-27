#' Calculate a relatedness matrix
#'
#' This function allows you to generate an n by n genetic relatedness matrix. If a numeric matrix is supplied, the RRM (Hayes, 2009) is used
#' and is computed XX'/p, where X is standardized. 
#' @param X A numeric matrix of genotypes (from *fully-imputed* data)
#' @param std Logical: should X be standardized? If you set this to FALSE, you should know exactly why that is appropriate... standardization is a best practice, and this will impact results.
#' @export
#' 
#' @examples 
#' RRM <- relatedness_mat(X = admix$X)
relatedness_mat <- function(X, std = TRUE){
  # NB: X is standardized as part of the RRM calculation
  if (std){
    rrm <- tcrossprod(ncvreg::std(X))/ncol(X)
  } else {
    rrm <- tcrossprod(X)/ncol(X)
  }
  
  return(rrm)
}

