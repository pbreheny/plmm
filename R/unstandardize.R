#' Unstandardize coefficient values
#'
#' This function allows you to unstandardize coefficient values based on attributes of X.
#' @param std_betas p x nlambda matrix of standardized coefficient path values.
#' @param X Original, non-standardized design matrix without an intercept column. Necessary for properly defining the dimensions of beta in cases where singular columns are present.
#' @param std_X Standardized design matrix. Must include the attributes 'center', 'scale', and 'nonsingular' corresponding to non-intercept variables.
#' @param intercept Logical flag for whether an intercept row is included in b, that is, a coefficient which accounts for the mean of the original y value on its rotated scale.
#' @export
#' 
#' @examples 
#' \dontrun{
#' admix$K <- relatedness_mat(admix$X) # create an estimated covariance matrix 
#' fit <- plmm(X = admix$X, y = admix$y, K = admix$K)
#' # TODO:need example here
#' }
#' 

unstandardize <- function(std_betas, X, std_X, intercept = TRUE) {
  # mark those columns which are non-singular
  # Note: singular columns will be assigned a beta value of 0
  ns <- attr(std_X, 'nonsingular')
  
  if (intercept){
  # partition the beta values into intercept and non-intercept parts
    a <- std_betas[1, , drop = FALSE] 
    b <- std_betas[-1, , drop=FALSE]
    
    # initialize beta with zeros; the + 1 is for the intercept
    un_std_beta <- matrix(0, nrow = (ncol(X) + 1), ncol = ncol(std_betas))
    
    # unscale the beta values for nonsingular columns
    unscaled_b <- b/attr(std_X, 'scale')[ns]
    
    # fill in the un-transformed values
    un_std_beta[ns+1,] <- unscaled_b # again, the + 1 is for the intercept
    un_std_beta[1,] <- a - crossprod(attr(std_X, 'center')[ns], unscaled_b)
    
  } else {
    un_std_beta <- matrix(0, nrow = ncol(X), ncol = ncol(b))
    unscaled_b <- b / scale
    un_std_beta[ns, ] <- unscaled_b
  }
  return(un_std_beta)
}
