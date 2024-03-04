#' helper function to implement MCP penalty
#'
#' @param z a vector representing the solution over active set at each feature
#' @param l1 upper bound (on beta)
#' @param l2 lower bound (on beta)
#' @param gamma The tuning parameter of the MCP penalty 
#' @param v the 'xtx' term 
#' @keywords internal
MCP <- function(z, l1, l2, gamma, v) {
  s <- rep(0, length(z))
  ret <- rep(NA_integer_, length(z))
  
  s[z > 0] <- 1
  s[z < 0] <- -1
  
  ret[abs(z) <= l1] <- 0
  ret[abs(z) <= gamma*l1*c(1 + l2)] <- s*(abs(z) - l1)/(v*c(1 + l2 - 1/gamma))
  ret[is.na(ret)] <- z/(v*(1 + l2))
  
  return(ret)
}

#' helper function to implement SCAD penalty
#'
#' @param z solution over active set at each feature
#' @param l1 upper bound
#' @param l2 lower bound
#' @param gamma The tuning parameter of the SCAD penalty 
#' @param v the 'xtx' term 
#'
#' @keywords internal
SCAD <- function(z, l1, l2, gamma, v) {
  s <- rep(0, length(z))
  ret <- rep(NA_integer_, length(z))

  s[z > 0] <- 1
  s[z < 0] <- -1
  
  ret[abs(z) <= l1] <- 0
  ret[abs(z) <= (l1 *c(1 + l2) + l1)] <- s * c(abs(z) - l1) / (v * c(1 + l2))
  ret[abs(z) <= gamma * l1 * c(1 + l2)] <- s * c(abs(z) - gamma * l1 / (gamma - 1)) / c(v * c(1 - 1 / (gamma - 1) + l2))
  
  ret[is.na(ret)] <- z / (v * c(1 + l2))
  
  return(ret)
}

#' helper function to implement lasso penalty
#'
#' @param z solution over active set at each feature
#' @param l1 upper bound
#' @param l2 lower bound 
#' @param v the 'xtx' term 
lasso <- function(z, l1, l2, v) {
  s <- rep(0, length(z))
  ret <- rep(NA_integer_, length(z))
  
  s[z > 0] <- 1
  s[z < 0] <- -1
  
  ret[abs(z) <= l1] <- 0
  ret[is.na(ret)] <- s * c(abs(z) - l1) / c(v * c(1 + l2))
  
  return(ret)
}



