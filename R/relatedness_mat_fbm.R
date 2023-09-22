#' A function to calculate the default realized relatedness matrix for a FBM
#'
#' @param X a FBM object as created by get_data()$X ()
#' @param ns a vector of values indicating the 
#' @param ... Other optional arguments (like ncores = ...)
#'
#' @return rrm The realized relatedness matrix, (1/p)*(XX^T)
#' @keywords internal
#'
relatedness_mat_fbm <- function(X, ns, ...){
  # cross product 
  cprod <- bigstatsr::big_tcrossprodSelf(X,
                                       # TODO: figure out why the line below makes the function return all NA values
                                       # fun.scaling = bigstatsr::as_scaling_fun(scale.col = rep(X$ncol, length(ns)),
                                       #                                         center.col = rep(0, length(ns))),
                                       ind.col = ns,
                                       ...)
  rrm <- FBM(cprod$nrow, cprod$ncol)
  # scale by p, the number of columns in the design matrix (including constant features)
  bigstatsr::big_apply(X = cprod,
                       a.FUN = function(X, ind, p, res){
                         res[,ind] <- sweep(x = X[,ind],
                                            MARGIN = 2,
                                            STATS = p,
                                            FUN = "/")},
                       a.combine = cbind,
                       p = X$ncol,
                       res = rrm)
  
  
  return(rrm)
}