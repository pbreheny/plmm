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
  
  ## setup 'a' (beta from previous iteration)----------------
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
    v <- bigstatsr::big_apply(X,
                   a.FUN = function(X, ind, n){
                     sum(X[,ind]^2)/n
                   },
                   a.combine = c,
                   ind = bigstatsr::cols_along(X),
                   ncores = bigstatsr::nb_cores(),
                   n = X$nrow)
    
  } else {
    v <- xtx 
  }
  
  
  # setup z 
  z <- rep(0, p)
  
  # take gaussian loss of y
  gl <- sum(y^2)
  # TODO: maybe move to C in future version
  # gl <- .Call("g_loss", y, n)  
  sd_y <- sqrt(gl/n)
  
  maxchange <- rep(0, length(y)) 
  
  # fit the model ------------------------------
  
  ## solve over the active set -----------------------------------------
  if(length(active) != 0){
    cat("\nSolving over active features")
    z[active] <- bigstatsr::big_apply(X = X,
                         a.FUN = function(X, ind, r, n, v, a){
                           # cp = cross product 
                           cp <- apply(X = X[,ind],
                                       MARGIN = 2,
                                       FUN = function(j){
                                         crossprod(j, r)/n
                                       })
                           
                           cp + (v[ind]*a[ind]) -> foo
                         },
                         a.combine = c,
                         ind = active, # this is the key line here
                         ncores = bigstatsr::nb_cores(),
                         r = r_new,
                         n = n,
                         v = drop(xtx),
                         a = a)
    
    # update beta ind, lam, multiplier, alpha, b, z, gamma, v, penalty
    b <- update_beta(ind = active, lam = lambda, alpha = alpha, b = b,
                     z = z, gamma = gamma, v = v, penalty = penalty,
                     multiplier = multiplier)
    
    # update r 
    r_new <- update_r(r = r_new, X = X, shift = b - a, ind = active)
    
    # update maxchange
    check1 <- abs(b - a)*sqrt(v) > maxchange
    maxchange[check1] <- (abs(b - a)*sqrt(v))[check1]
    
  } 
  
  # update coefficient estimates 
  a[active] <- b[active]
  # check for convergence 
  if(any(maxchange > eps*sd_y)){
  warning("\nConvergence issue")
    
  }
  
  ## scan for violations -----------------------------------------------
  inactive <- which(a == 0)
  violations <- 0
  z[inactive] <- bigstatsr::big_apply(X = X,
                       a.FUN = function(X, ind, r, n){
                         cp <- crossprod(X[,ind], r) # note: no 'v' and 'a' here
                         cp/n 
                       },
                       a.combine = c,
                       ind = inactive, # note: this is the key change from above
                       ncores = bigstatsr::nb_cores(),
                       r = r_new,
                       n = n)
  
  
  # update beta 
  cat("\nSolving over inactive features")
  b <- update_beta(ind = inactive, lam = lambda, alpha = alpha, b = b,
                   z = z, gamma = gamma, v = v, penalty = penalty,
                   multiplier = multiplier)
  
  # if something enters, update active set & residuals
  nz_beta <- which(b != 0)
  active[nz_beta] <- 1
  r_new <- update_r(r_new, X, shift = b, ind = inactive)
  
  # update coefficient estimates
  a[inactive] <- b[inactive]
  
  # update violation count 
  violations <- violations + length(nz_beta) 
  
  
# return result to plmm_fit
res <- list(
  beta = b,
  loss = sum(r^2),
  resid = r
)

return(res)



}

