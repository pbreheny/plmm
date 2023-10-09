#' big_std
#'A function to standardize a matrix stored as an FBM object 
#'
#' @param X An file-backed design matrix 
#' @param center A vector of centering values
#' @param scale A vctor of scaling values 
#' @param ns A vector with the indicies marking the Non-Singular columns of X
#' @param fbm Logical: should the standardized matrix returned be an FBM object? Defaults to TRUE. 
#' @param ... Other arguments to `bigstatsr::big_apply()`
#'
#' @return A standardized matrix of class FBM
#' @keywords internal
big_std <- function(X, center = NULL, scale, ns, fbm = TRUE){
  if(is.null(center)){
    centered_X <- X
  } else {
    # allocate space 
    centered_X <- bigstatsr::FBM(X$nrow, X$ncol)
    bigstatsr::big_apply(X = X,
                         a.FUN = function(X, ind, center, res){
                           res[,ind] <- sweep(x = X[,ind],
                                              MARGIN = 2,
                                              STATS = center[ind],
                                              FUN = "-")},
                         a.combine = cbind,
                         ncores = bigstatsr::nb_cores(),
                         center = center,
                         res = centered_X)
  }
  
  
  scaled_X <- bigstatsr::FBM(centered_X$nrow, centered_X$ncol)
  bigstatsr::big_apply(X = centered_X, 
                       a.FUN = function(X, ind, scale, res){
                         res[,ind] <- sweep(x = X[,ind],
                               MARGIN = 2,
                               STATS = scale[ind],
                               FUN = "/")},
                       a.combine = cbind,
                       # NB: only scale the nonsingular columns
                       ind = ns, 
                       ncores = bigstatsr::nb_cores(),
                       scale = scale,
                       res = scaled_X)
  
  
  
  return(scaled_X)
  
}
