#' helper function to implement MCP penalty
#'
#' @param z solution over active set at each feature
#' @param l1 upper bound (on beta)
#' @param l2 lower bound (on beta)
#' @param gamma The tuning parameter of the MCP penalty 
#' @param v the 'xtx' term 
#' @keywords internal
MCP <- function(z, l1, l2, gamma, v) {
  s <- 0
  if (z > 0) s <- 1
  else if (z < 0) s <- -1
  if (abs(z) <= l1) return(0)
  else if (abs(z) <= gamma * l1 * (1 + l2)) return(s * (abs(z) - l1) / (v * (1 + l2 - 1 / gamma)))
  else return(z / (v * (1 + l2)))
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
  s <- 0
  if (z > 0) s <- 1
  else if (z < 0) s <- -1
  if (abs(z) <= l1) return(0)
  else if (abs(z) <= (l1 * (1 + l2) + l1)) return(s * (abs(z) - l1) / (v * (1 + l2)))
  else if (abs(z) <= gamma * l1 * (1 + l2)) return(s * (abs(z) - gamma * l1 / (gamma - 1)) / (v * (1 - 1 / (gamma - 1) + l2)))
  else return(z / (v * (1 + l2)))
}

#' helper function to implement lasso penalty
#'
#' @param z solution over active set at each feature
#' @param l1 upper bound
#' @param l2 lower bound 
#' @param v the 'xtx' term 
lasso <- function(z, l1, l2, v) {
  s <- 0
  if (z > 0) s <- 1
  else if (z < 0) s <- -1
  if (abs(z) <= l1) return(0)
  else return(s * (abs(z) - l1) / (v * (1 + l2)))
}



