#' a function to update beta
#'
#' @param active Indices of active features 
#' @param alpha 
#' @param b 
#' @param z 
#' @param gamma 
#' @param v 
#' @param lam A single value of lambda 
update_beta <- function(active, lam, alpha, b, z, gamma, v, penalty){
  browser()
  for (j in 1:active){
    l1 = lam * multiplier[j] * alpha
    l2 = lam * multiplier[j] * (1-alpha)
    
    if (penalty == "MCP"){b[j] <- MCP(z[j], l1, l2, gamma, v[j])} 
    if (penalty == "SCAD") {b[j] <- SCAD(z[j], l1, l2, gamma, v[j])}
    if (penalty == "lasso") {b[j] <- lasso(z[j], l1, l2, v[j])}
  }
}

#' a function to update r 
#'
#' @param n Number of observations 
#' @param active Indices of active features 
#' @param b Vector of estimated coefficients
#' @param a Vector of initial values
#' @param r Vector of residuals
update_r <- function(n, active, b, a, r){
  for(j in 1:active){ 
    shift <- b[j] - a[j]
    
    if (shift !=0) {
      for (i in 1:n){
        r[i] = r[i] - shift*X[j*n+i]
      }
      
      if(abs(shift)*sqrt(v[j]) > maxChange)maxChange <- abs(shift)*sqrt(v[j])
    }
    
  }
  
}
