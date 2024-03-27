#' Predict method for plmm class
#'
#' @param object An object of class \code{plmm}.
#' @param newX Design matrix used for computing predicted values if requested. 
#'  **Columns must be named!**
#' @param type A character argument indicating what type of prediction should be 
#' returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param lambda A numeric vector of regularization parameter \code{lambda} values 
#' at which predictions are requested.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which 
#' predictions are required. By default, all indices are returned.
#' @param X Original design matrix (not including intercept column) 
#' from object. Required if \code{type == 'blup'} and object is too large to be 
#' returned in `plmm` object. Again, **columns must be named!**
#' @param y Original continuous outcome vector from object. 
#' Required if \code{type == 'blup'}. 
#' @param K An optional list or matrix as returned by `choose_K()`. 
#' @param ... Additional optional arguments
#' 
#' @details
#' Define beta-hat as the coefficients estimated at the value of lambda that minimizes cross-validation error (CVE). Then options for `type` are as follows: 
#' 
#'  * 'response' (default): uses the product of newX and beta-hat to predict new values of the outcome. This does not incorporate the correlation structure of the data. 
#'  For the stats folks out there, this is simply the linear predictor. 
#'  
#'  * 'blup' (acronym for Best Linear Unbiased Predictor): adds to the 'response' a value that represents the esetimated random effect. This addition is a way of incorporating 
#'  the estimated correlation structure of data into our prediction of the outcome. 
#'  
#'  * 'coefficients': returns the estimated beta-hat 
#'  
#'  * 'vars': returns the _indicies_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept. 
#'  
#'  * 'nvars': returns the _number_ of variables (e.g., SNPs) with nonzero coefficients at each value of lambda. EXCLUDES intercept. 
#' 
#' 
#' @rdname predict.plmm
#' @export
#'
#' @examples 
#' \dontrun{
#' # fit a model 
#' set.seed(123)
#' train_idx <- sample(1:nrow(admix$X), 100) # shuffling is important here! Keeps test and train groups comparable.
#' train <- list(X = admix$X[train_idx,], y = admix$y[train_idx])
#' test <- list(X = admix$X[-train_idx,], y = admix$y[-train_idx])
#' fit <- plmm(X = train$X, y = train$y, K = relatedness_mat(train$X))
#' 
#'  # make predictions for all lambda values 
#'  pred1 <- predict(object = fit, newX = test$X, type = "blup", X = train$X, y = train$y)
#'
#'  # make predictions for a select number of lambda values
#'  # use cv to choose a best lambda
#'  cvfit <- cv.plmm(X = train$X, y = train$y) 
#'  pred2 <- predict(object = fit, newX = test$X, type = "blup", idx=cvfit$min,
#'   X = train$X, y = train$y)
#'  pred3 <- predict(object = fit, newX = test$X, type = "lp", idx=cvfit$min) 
#'  
#'   # examine prediction error 
#'   summary(crossprod(test$y - pred2))  # for BLUP method 
#'   summary(crossprod(test$y - pred3)) # for LP method 
#'   
#'   # which was best prediction?
#'   mspe <- apply(pred1, 2, function(c){crossprod(test$y - c)})
#'   which.min(mspe) # not anywhere near the 'best lambda' chosen in cross validation...
#'   mspe[which.min(mspe)]
#'   min(cvfit$cve)
#'  }
#'  
#'  
predict.plmm <- function(object,
                         newX,
                         type=c("lp", "coefficients", "vars", "nvars", "blup"),
                         lambda,
                         idx=1:length(object$lambda),
                         X,
                         y,
                         K = NULL,
                         ...) {
  
  type <- match.arg(type)
  beta_vals <- coef.plmm(object, lambda, which=idx, drop=FALSE) # includes intercept 
  p <- object$p 
  n <- object$n 

  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)

  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept
  
  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  
  a <- beta_vals[1,]
  b <- beta_vals[-1,,drop=FALSE]

  Xb <- sweep(newX %*% b, 2, a, "+")
  # Xb <- cbind(1, newX) %*% beta_vals # old way 
  
  if (type=="lp") return(drop(Xb))
  
  if (type == "blup"){ # assuming eta of X and newX are the same 
    
    if (missing(X)) stop("The design matrix is required for BLUP calculation. Please supply the no-intercept design matrix to the X argument.") 
    if (missing(y) & is.null(object$y)) stop("The vector of outcomes is required for BLUP calculation. Please either supply it to the y argument, or set returnX=TRUE in the plmm function.")
    
    if (!is.null(object$y)) y <- object$y
    
    if (is.null(K)){
      V11 <- object$eta*relatedness_mat(X, std = FALSE) + (1 - object$eta)*diag(object$n)
    } else {
      V11 <- object$eta*K + (1 - object$eta)*diag(object$n)
    }

    # standardize (this won't affect predicted y)
    # std_X <- ncvreg::std(X)
    # std_newX <- ncvreg::std(newX)
    # std_X_details <- list(ns = attr(std_X, 'nonsingular'),
    #                       scale = attr(std_X, 'scale'),
    #                       center = attr(std_X, 'center'))
    
    # check dimensions -- must have same number of features
    # even from the same dataset, this can be an issue with features 'becoming' constant...
    p1 <- ncol(X)+1
    p2 <- ncol(newX)+1
    p3 <- nrow(beta_vals)
    if (!identical(p1, p2) | !identical(p1, p3)) {
      shared_p <- intersect(colnames(X), colnames(newX)) # this is why we need column names
      X <- X[,shared_p]
      newX <- newX[,shared_p]
      beta_vals <- beta_vals[c("(Intercept)",shared_p),]
    }
    
    # cannot use U, S when computing V21, because nv=0. V is needed to restore X 
    V21 <- object$eta * (1/p) * tcrossprod(newX,X) # same as V21_check 
    Xb_old <- sweep(X %*% b, 2, a, "+")
    resid_old <- drop(y) - Xb_old
    
    ranef <- V21 %*% (chol2inv(chol(V11)) %*% resid_old)
    blup <- drop(Xb + ranef)


    # #################################
    # # used to check BLUP calculation #
    # #################################
    # aggregate X and newX to compute V
    # X_all <- rbind(std_X, std_newX)
    # c(S_all, U_all) %<-% svd(X_all, nv = 0) # D, U
    # S_all <- S_all^2 / p

    # assuming newX has the same eta as X
    # eta_all <- object$eta
    # Vhat_all <- eta_all * tcrossprod(U_all %*% diag(S_all), U_all) + (1-eta_all)*diag(nrow(U_all))
    # V21_check <- Vhat_all[-c(1:n), 1:n, drop = FALSE]
    # V11_check <- Vhat_all[1:n, 1:n, drop = FALSE]

    # ranef_check <- V21_check %*% chol2inv(chol(V11_check)) %*% (drop(y) - cbind(1, std_X) %*% beta_vals)
    # print(eta)

    # blup_check <- Xb + ranef_check
    
    # frob_distance <- (Matrix::norm(blup_check) - Matrix::norm(blup))/Matrix::norm(blup_check)
    # if (abs(frob_distance) > 0.01) stop("\nThe two calculations of the BLUP do not match")
    

    return(blup)
  }

}
