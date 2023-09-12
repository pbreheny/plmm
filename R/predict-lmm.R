#' Predict method for lmm class
#'
#' @param object An object of class \code{plmm}.
#' @param newX Design matrix used for computing predicted values if requested.
#' @param type A character argument indicating what type of prediction should be returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param X Optional argument. Original design matrix (not including intercept column) from object. Required if \code{type == 'blup'} and object is too large to be returned in `plmm` object.
#' @param y Optional argument. Original continuous outcome vector from object. Required if \code{type == 'blup'}.
#' @param ... Additional optional arguments
#' 
#' @details
#' Define beta-hat as the coefficients estimated at the value of lambda that minimizes cross-validation error (CVE). Then options for `type` are as follows: 
#' 
#'  * 'lp' (default): uses the product of newX and beta-hat to predict new values of the outcome. This does not incorporate the correlation structure of the data. 
#'  For the stats folks out there, this is simply the linear predictor. 
#'  
#'  * 'blup' (acronym for Best Linear Unbiased Predictor): adds to the 'response' a value that represents the esetimated random effect. This addition is a way of incorporating 
#'  the estimated correlation structure of data into our prediction of the outcome. 
#'  
#'  * 'coefficients': returns the estimated beta-hat 
#' 
#' @rdname predict.lmm
#' @export
#'
#' @examples 
#' \dontrun{
#' # fit a model 
#' fit <- lmm(X = pedigree$X, y = pedigree$clinical$y)
#' predict.lmm(fit, type = "coefficients")
#' predict.lmm(fit, type = "lp", newX = pedigree$X)
#' predict.lmm(fit, type = "blup", X = pedigree$X, newX = pedigree$X, y = pedigree$clinical$y)
#' }
#' 
predict.lmm <- function(object,
                         newX,
                         type=c("lp", "coefficients", "blup"),
                         X,
                         y,
                         ...) {
  
  type <- match.arg(type)
  beta_vals <- coef.lmm(object) # includes intercept 
  p <- object$ncol_X 
  n <- object$nrow_X 
  
  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)
  
  if(type %in% c('blup', 'lp') & missing(newX)){
    stop("Predictions requested but no new data supplied; please specify newX argument.")
  }
  
  Xbeta <- cbind(1, newX) %*% beta_vals
  if (type=="lp") return(drop(Xbeta))
  
  if (type == "blup"){ # assuming eta of X and newX are the same 
    
    message("The BLUP option is a newer development. If you notice anything unusal or encounter any problems, please file an issue on our GitHub repo.") 
    
    if (missing(X)) stop("The design matrix is required for BLUP calculation. Please supply the no-intercept design matrix to the X argument.") 
    if (missing(y) & is.null(object$y)) stop("The vector of outcomes is required for BLUP calculation. Please either supply it to the y argument, or set returnX=TRUE in the plmm function.")
    # if (!(c("S", "U") %in% names(object))) stop("SVD results are required for BLUP calculation. Use 'svd_details = TRUE' in the 'plmm' function.")
    
    if (!is.null(object$y)) y <- object$y
    
    ###################
    # used to compute #
    ################### 
    
    V11 <- object$Vhat # same as V11_check 
    
    # cannot use U, S when computing V21, because nv=0. V is needed to restore X 
    V21 <- object$eta * (1/p) * tcrossprod(newX, X) # same as V21_check 
    
    ranef <- V21 %*% chol2inv(chol(V11)) %*% (drop(y) - cbind(1, X) %*% beta_vals)
    blup <- drop(Xbeta + ranef)
    
    # #################
    # # used to check #
    # #################
    # 
    # # aggregate X and newX to compute V 
    # X_all <- rbind(X, newX)
    # c(S_all, U_all) %<-% svd(X_all, nv = 0) # D, U
    # S_all <- S_all^2 / p
    # 
    # # assuming newX has the same eta as X 
    # eta_all <- object$eta
    # Vhat_all <- eta_all * tcrossprod(U_all %*% diag(S_all), U_all) + (1-eta_all)*diag(nrow(U_all)) 
    # V21_check <- Vhat_all[-c(1:n), 1:n, drop = FALSE] 
    # V11_check <- Vhat_all[1:n, 1:n, drop = FALSE] 
    # # V11 <- object$estimated_V
    # 
    # ranef_check <- V21_check %*% chol2inv(chol(V11_check)) %*% (drop(y) - cbind(1, X) %*% beta_vals)
    # # print(eta) 
    # 
    # blup_check <- Xbeta + ranef_check
    
    return(blup)
  }
  
}




