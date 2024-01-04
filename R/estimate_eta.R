#' Estimate eta (to be used in rotating the data)
#' This function is called internally by \code{plmm()}
#' @param s The singular values of K, the realized relationship matrix
#' @param U The left-singular vectors of the *standardized* design matrix
#' @param y Continuous outcome vector.
#' @export
#' @keywords internal
estimate_eta <- function(s, U, y){

  # optimize over \hat\beta_0 and \eta
  opt <- optim(par = c(0.1, 0.01),
                  fn = null_model_nll,
                  y = y,
                  s = s,
                  U = U,
                  method = "BFGS")
  
  return(list(eta = c(opt$par[1]),
              beta0 = opt$par[2]))
}

#' a helper function for 2-dim optimization
#' @param params A vector of 2 elements, corresponding to eta and beta0 
#' The latter is the coefficient of the null model
#' @param y A vector of outcomes 
#' @param U The left singular vectors of data X 
#' @param s the vector of singular values of data X 
#' 
#' @keywords internal
null_model_nll <- function(params, y, U, s){
  
  # name parameters
  eta <- params[1]
  beta0 <- params[2]
  
  
  # create intercept (for null model, this is the only predictor)
  intcpt <- rep(1, length(y))
  
  # get \Sigma^2_{-1/2} piece
  w <- (eta * s + (1 - eta))^(-1/2)
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  rot_y <- wUt %*% y 
  rot_intcpt <- wUt %*% intcpt
  # browser()
  # wUt_svd <- tcrossprod(wUt) |> svd()
  # distribution of null model 
  res <- mvtnorm::dmvnorm(x = drop(rot_y),
               mean = drop(rot_intcpt*beta0),
               sigma = tcrossprod(wUt),
               log = TRUE)
  
  ret <- sum(res)
  
  return(-1*ret) # want to use minimization in optim()

}
