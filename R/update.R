#' a function to update beta
#'
#' @param ind     Indices of features (e.,g., active features)
#' @param alpha   Tuning paramter controlling the ridge component of penalty, as 
#'                in `ncvreg()`; default is 1 (meaning no ridge penalty)
#' @param b       Beta coefficient values 
#' @param z       The solution over the active set 
#' @param gamma   Tuning parameter of the MCP/SCAD penalty, as in `ncvreg()`; 
#'                default is 3 for MCP and 3.7 for SCAD.
#' @param v       
#' @param lam     A single value of lambda 
update_beta <- function(ind, lam, multiplier, alpha, b, z, gamma, v, penalty){
  
  l1 <- (lam*alpha)*multiplier[ind]
  l2 <- (lam*(1-alpha))*multiplier[ind]
  
  if(any(is.na(l1))){
    stop("\nMissing values in l1 arg of update beta")
  }
  
  if (penalty == "MCP"){b[ind] <- MCP(z[ind], l1, l2, gamma, v[ind])} 
  if (penalty == "SCAD") {b[ind] <- SCAD(z[ind], l1, l2, gamma, v[ind])}
  if (penalty == "lasso") {b[ind] <- lasso(z[ind], l1, l2, v[ind])}
  
  return(b)
}

#' a function to update r 
#' @param r Vector of residuals
#' @param ind Indicies for which features to use in update (e.g., indicies for active set)
#' @param X A file-backed data frame 
#' @param shift A vector of the same length of `r`, indicating how much to shift/change `X`
update_r <- function(r, X, shift, ind){
  change <- rep(0, length(r))
  X <- bigstatsr::big_copy(X, ind.col = ind)
  change <- bigstatsr::big_prodVec(X, shift[ind])
  r - change
  # return(r)
}
  

