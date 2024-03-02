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
  browser()
    l1 <- lam * multiplier[active] * alpha
    l2 <- lam * multiplier[active] * (1-alpha)
    
    if (penalty == "MCP"){b[ind] <- MCP(z[ind], l1, l2, gamma, v[ind])} 
    if (penalty == "SCAD") {b[ind] <- SCAD(z[ind], l1, l2, gamma, v[ind])}
    if (penalty == "lasso") {b[ind] <- lasso(z[ind], l1, l2, v[ind])}
}

#' a function to update r 
#'
#' @param n Number of observations 
#' @param ind Indices of features 
#' @param b Vector of estimated coefficients
#' @param a Vector of initial values
#' @param r Vector of residuals
#' @param X A file-backed data frame 
update_r <- function(n, ind, b, a, r, X){
  shift <- b[ind] - a[ind]
  
  if (shift !=0) {
    for (i in 1:n){
      r[i] = r[i] - shift*X[j*n+i]
    }
    
    if(abs(shift)*sqrt(v[j]) > maxChange)maxChange <- abs(shift)*sqrt(v[j])
  }
  
}
  

