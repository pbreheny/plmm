#' std_fbm
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
big_std <- function(X, center, scale, ns, fbm = TRUE){
  centeredX <- bigstatsr::big_apply(X = X,
                                    a.FUN = function(X, ind, center){
                                      sweep(x = X[,ind],
                                            MARGIN = 2,
                                            STATS = center[ind],
                                            FUN = "-")},
                                    a.combine = cbind,
                                    ncores = bigstatsr::nb_cores(),
                                    center = center)
  centeredXfbm <- bigstatsr::as_FBM(centeredX)

  scaledX <- bigstatsr::big_apply(X = centeredXfbm, 
                                  a.FUN = function(X, ind, scale){
                                    sweep(x = X[,ind],
                                          MARGIN = 2,
                                          STATS = scale[ind],
                                          FUN = "/")},
                                  a.combine = cbind,
                                  # NB: only scale the nonsingular columns
                                  ind = ns, 
                                  ncores = bigstatsr::nb_cores(),
                                  scale = scale)
  if(fbm){
    return(bigstatsr::as_FBM(scaledX))
  } else {
    return(scaledX)
}
}
