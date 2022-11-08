#' Compute transformed X and y for subsequent plmm rescaling and fitting
#'
#' This function allows you to fit a linear mixed model via non-convex penalized maximum likelihood.
#' @param X Design matrix. May include clinical covariates and other non-SNP data.
#' @param y Continuous outcome vector.
#' @param K Estimated or known similarity matrix. By default, K is the realized relationship matrix, \eqn{\frac{1}{p}XX^T}, where \eqn{p} is the number of columns in X
#' @param eta_star Optional arg to input a specific eta term rather than estimate it from the data. If v is a known matrix, this should be 1.
#' @param k An integer between 1 and \code{nrow(K)} indicating the number of singular values requested *if* package \code{RSpectra} is installed. Defaults to NULL. 
#' @importFrom zeallot %<-%
#' @export
#' 
#' @examples 
#' std_X <- scale(admix$X)
#' K <- relatedness_mat(std_X)
#' rotated_dat <- rotate_data(std_X, admix$y, K)
#' rotated_dat2 <- rotate_data(std_X, admix$y, K = K, k = 5)

rotate_data <- function(X, y, K = NULL, eta_star, k = NULL){
  # Coersion
  S <- U <- UU <- eta <- NULL
  
  # check to see if RSpectra is installed 
  RSpectra <- ifelse(system.file(package = "RSpectra") == "",
                     FALSE, 
                     TRUE)
  

  # Calculate RRM
  if(is.null(K)){
    p <- ncol(X)
    K <- relatedness_mat(X)
  }
  
  ## Calculate U, S, eta
  # case 1: estimate eta if needed
  if (missing(eta_star)){
    if(RSpectra & !is.null(k)){
      c(S, U, eta) %<-% plmm_null(K = K, y = y, k = k) 
    } else {
      c(S, U, eta) %<-% plmm_null(K = K, y = y)
    }
    
    # case 2: eta is supplied as an argument 
  } else {
    if(RSpectra & !is.null(k)){
      c(S, U) %<-% Rspectra::svds(K, nv = 0, k = k)[1:2]
    } else {
      c(S, U) %<-% svd(K)[1:2]
    }
    
      # NB: nv=0 avoids the calculation of right singular vectors, which we will not use
      eta <- eta_star
    }
  # Construct W 
  W <- diag((eta * S + (1 - eta))^(-1/2))
  

  ## Rotate data
  SUX <- W %*% crossprod(U, cbind(1, X)) # add column of 1s for intercept
  SUy <- drop(W %*% crossprod(U, y))
  val <- list(SUX = SUX,
              SUy = SUy,
              eta = eta,
              U = U,
              S = S)
  return(val)
}
