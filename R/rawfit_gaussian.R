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
  
  # error checking ---------------------------
  # TODO: come back and add steps similar to ncvreg::ncvfit()
  
  # setup steps -------------------------------
  n <- length(y)
  p <- X$ncol # TODO: check this against ncvreg::rawfit_gaussian.c
  
  ## initialize 'b' (new beta) --------------------
  b <- rep(0, p)
  
  ## setup 'a' (use beta from previous iteration)----------------
  a <- init
  
  # initialize indicators for active set 
  active <- which(a != 0) 
  
  ## setup r (the residuals) ----
  r_new <- as.double(r)
  # address the case below where r could be NA
  if (is.na(r_new[1])) {
    bigstatsr::big_apply(X = X,
                         a.FUN = function(X, ind, y, a, res){
                           res[ind] <- y - crossprod(X[ind,], a)
                         },
                         a.combine = c,
                         ind = bigstatsr::rows_along(X),
                         ncores = bigstatsr::nb_cores(),
                         a = a,
                         res = r_new,
                         y = y)
  } 
  
  ## setup v (sums of squares) ----
  v <- rep(NA_integer_, )
  if (is.na(xtx[1])){
    # create v for case where xtx is NA 
    # GOAL: function(X, n, j) {
    #   nn <- n * j
    #   val <- sum(X[nn + (1:n)]^2)
    #   return(val)
    # }
    # TODO: translate the above into a big_apply() call
    #  big_apply()
    
  } else {
    v <- xtx 
  }
  
  
  # setup z 
  z <- rep(NA_integer_, p)
  
  # take gaussian loss of y
  gl <- sum(y^2)
  # gl <- .Call("g_loss", y, n) # TODO: maybe move to C in future version 
  sd_y <- sqrt(gl/n)
  maxchange <- 0 
  
  # fit the model ------------------------------
  
  ## solve over the active set -----------------------------------------
  if(length(active) != 0){
    bigstatsr::big_apply(X = X,
                         a.FUN = function(X, ind, r, n, v, a, res){
                           cp <- crossprod(X[,ind], r)
                           res[ind] <- cp/n + (v[ind]*a[ind])
                         },
                         a.combine = c,
                         ind = active, # this is the key line here
                         ncores = bigstatsr::nb_cores(),
                         r = r_new,
                         n = n,
                         v = drop(xtx),
                         a = a,
                         res = z)
    
    # update beta 
    b <- update_beta(active, lambda, alpha, b, z, gamma, v, penalty)
    
    # update r 
    r <- update_r(n, active, b, a, r)
    
    
  } 
  
  # check for convergence 
  for (j in 1:p){
    a[j] <- b[j]
    if(maxchange < eps*sd_y) break;
  }
  
  ## scan for violations -----------------------------------------------
  nonactive <- which(a == 0)
  violations <- 0
  z <- bigstatsr::big_apply(X = X,
                       a.FUN = function(X, ind, r, n, res){
                         cp <- crossprod(X[,ind], r) # note: no 'v' and 'a' here
                         res[ind] <- cp/n 
                       },
                       a.combine = c,
                       ind = nonactive, # note: this is the key change from above
                       ncores = bigstatsr::nb_cores(),
                       r = r_new,
                       n = n,
                       res = z)
  
  
  # update beta 
  b <- update_beta(nonactive, lam = lambda, alpha, b, z, gamma, v)
  browser()
  # if something enters, update active set & residuals
  nz_beta <- which(b != 0)
  active[nz_beta] <- 1
  bigstatsr::big_apply(X = X,
                       a.FUN = function(X, ind, b, res){
                         res[ind] <- res[ind] - b[ind]*X[ind,]},
                       a.combine = c,
                       ind = rows.along(X),
                       b = b,
                       res = r)
  
  
# return result to plmm_fit
res <- list(
  beta = b,
  loss = sum(r^2),
  resid = r
)

return(res)



}

