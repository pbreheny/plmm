#' a function to fit a gaussian model for a single value of lambda on a filebacked X 
#' 
#' @param X 
#' @param y
#' @param init
#' @param r
#' @param xtx
#' @param penalty
#' @param lambda
#' @param eps
#' @param max_iter
#' @param gamma
#' @param multiplier
#' @param alpha
#' 
#' @keywords internal
#' 
rawfit_gaussian <- function(X, y, init, r, xtx, penalty, lambda, eps, max_iter, 
                            gamma, multiplier, alpha){
  
  browser() # PICK UP HERE!
  
  # error checking ---------------------------
  # TODO: come back and add steps similar to ncvreg::ncvfit()
  
  # setup steps -------------------------------
  n <- length(y)
  p <- X$ncol # TODO: check this against ncvreg::rawfit_gaussian.c
  
  # initialize 'b' (new beta)
  b <- rep(0, p)
  
  # setup 'a' (use beta from previous iteration)
  a <- init
  
  # initialize indicators for active set 
  active <- which(a != 0) 
  
  # setup r (the residuals)
  r_new <- r
  # TODO: address the case below where r could be NA
  # if(is.na(r[1])){
  #   r <- y
  #   
  # }
  
  
  
  # setup v (sums of squares)
  # TODO: come back and create v for case where xtx is NA 
  
  # setup z 
  z <- rep(NA_integer_, p)
  
  # take gaussian loss of y
  gl <- .Call("g_loss", y, n)
  sd_y <- sqrt(gl/n)
  maxchange <- 0 
  # fit the model ------------------------------
  
  # solve over the active set 
  bigstatsr::big_apply(X = X,
                       a.FUN = function(X, ind, r, n, v, a, res){
                         cp <- crossprod(X[,ind], r)
                         res[ind] <- cp/n + (v[ind]*a[ind])
                       },
                       a.combine = c,
                       ind = active,
                       ncores = bigstatsr::nb_cores(),
                       r = r_new,
                       n = n,
                       v = drop(xtx),
                       a = a,
                       res = z)
  
  # update beta 
  b <- update_beta(active, lam, alpha, b, z, gamma, v)

  # update r 
  r <- update_r(n, active, b, a, r)
  
  # check for convergence 
  for (j in 1:p){
    a[j] <- b[j]
    if(maxchange < eps*sdy) break;
  }
  
  # scan for violations -----
  nonactive <- which(a == 0)
  violations <- 0
  
  # solve over the inactive set 
  bigstatsr::big_apply(X = X,
                       a.FUN = function(X, ind, r, n, v, a, res){
                         cp <- crossprod(X[,ind], r)
                         res[ind] <- cp/n + (v[ind]*a[ind])
                       },
                       a.combine = c,
                       ind = nonactive, # note: this is the key change from above
                       ncores = bigstatsr::nb_cores(),
                       r = r_new,
                       n = n,
                       v = drop(xtx),
                       a = a,
                       res = z)
  
  
  # update beta 
  b <- update_beta(nonactive, lam, alpha, b, z, gamma, v)
  
  # if something enters, update active set & residuals
  nz_beta <- which(b != 0)
  active[nz_beta] <- 1
  # TODO: pick up here
  bigstatsr::big_apply(X = X,
                       a.FUN = function(X, ind, b, res){
                         res[ind] <- res[ind] - b[ind]*X[ind,]},
                       a.combine = c,
                       ind = rows.along(X),
                       b = b,
                       res = r)
  
  
  # return result to plmm_fit
  
  
  
  
  
}

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

#' a function to update beta
#'
#' @param active Indices of active features 
#' @param alpha 
#' @param b 
#' @param z 
#' @param gamma 
#' @param v 
#' @param lam A single value of lambda 
update_beta <- function(active, lam, alpha, b, z, gamma, v){
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


