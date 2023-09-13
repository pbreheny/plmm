#' Predict method for plmm class
#'
#' @param object An object of class \code{plmm}.
#' @param newX Design matrix used for computing predicted values if requested.
#' @param type A character argument indicating what type of prediction should be returned. Options are "lp," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param lambda A numeric vector of regularization parameter \code{lambda} values at which predictions are requested.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
#' @param X Optional argument. Original design matrix (not including intercept column) from object. Required if \code{type == 'blup'} and object is too large to be returned in `plmm` object.
#' @param y Optional argument. Original continuous outcome vector from object. Required if \code{type == 'blup'}.
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
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' 
#' # simulate new data 
#' newX <- sim_ps_x(n = nrow(admix$X), nJ = 4, p = ncol(admix$X),
#'  structureX = "independent", inbr = "heterogeneous", standardizeX = FALSE)
#'  
#'  # make predictions for all lambda values 
#'  pred1 <- predict(object = fit, newX = newX, type = "lp")
#'  
#'  # make predictions for a select number of lambda values 
#'  pred2 <- predict(object = fit, newX = newX, type = "lp", idx=98)
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
  p <- object$ncol_X 
  n <- object$nrow_X 
  
  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)

  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept
  
  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  
  Xbeta <- cbind(1, newX) %*% beta_vals
  if (type=="lp") return(drop(Xbeta))
  
  if (type == "blup"){ # assuming eta of X and newX are the same 
    
    if (missing(X)) stop("The design matrix is required for BLUP calculation. Please supply the no-intercept design matrix to the X argument.") 
    if (missing(y) & is.null(object$y)) stop("The vector of outcomes is required for BLUP calculation. Please either supply it to the y argument, or set returnX=TRUE in the plmm function.")
    
    if (!is.null(object$y)) y <- object$y
    
    ###################
    # used to compute #
    ################### 
    
    V11 <- v_hat(object, K)
    
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




# Appendix: Ideas from Anna's work 


# fit <- plmm(X = admix$X, y = admix$y, lambda = c(0.05, 0.01))
# head(predict(object = fit, newX = admix$X, type = 'lp', lambda = 0.05))
# head(predict(object = fit, newX = admix$X, type = 'vars'))
# predict(object = fit, newX = admix$X, type = 'nvars')
# \dontrun{
# predict.plmm(object = fit, newX = admix$X, type = 'blup', no_int_X = admix$X,
# y = admix$y, U = fit$U, S = fit$S, eta = fit$eta, covariance = crossprod(admix$X))
# }
# \dontrun{
# # Out of sample predictions are NOT improved by BLUP when X cannot be used to
# # estimate V
# n <- 210
# p <- 1000
# p1 <- 50
# SNR <- 10
#
# # Generate data (X unrelated to V)
# V <- Matrix::bdiag(matrix(0.5, n/3, n/3),
# matrix(0.5, n/3, n/3),
# matrix(0.5, n/3, n/3))
# V <- as.matrix(V)
# diag(V) <- 1
# set.seed(7)
# X1 <- matrix(rnorm(n * p, 0, 1), n, p)
# j <- 1:p
# s <- rep(1, times = p)
# b <- j <= p1
# b <- b * s
# covX <- stats::var(X1) * (n - 1)/n
# beta <- b * sqrt(SNR) / sqrt(drop(crossprod(b, covX) %*% b))
# Xbeta <- X1 %*% beta
# e <- MASS::mvrnorm(n=1, rep(0, n), V) # errors are structured
# y1 <- Xbeta + e
#
# # fit the model
# object <- cv.plmm(X1,
#                   y1,
#                   V,
#                   eta_star = 1,
#                   type = 'blup',
#                   penalty = "lasso",
#                   alpha = 1,
#                   standardizeX = TRUE,
#                   standardizeRtX = TRUE,
#                   rotation = TRUE,
#                   returnX = TRUE,
#                   nfolds = 5)
# # NOTE: in-sample predictions as above in CV are fine, because V supplies the
# # known covariance between the outcomes used in different folds.
#
# # generate an X2 that is similar to X1
# X2 <- X1 + rnorm(nrow(X1), 0, 0.0001)
# covX <- stats::var(X2) * (n - 1)/n
# beta <- b * sqrt(SNR) / sqrt(drop(crossprod(b, covX) %*% b))
# Xbeta <- X2 %*% beta
# e <- MASS::mvrnorm(n=1, rep(0, n), V)
# y2 <- Xbeta + e
#
# # linear predictor works well
# linear_predictor <- predict.plmm(object$fit,
# newX = X2,
# type = 'lp',
# lambda = object$lambda.min)
# drop(crossprod(linear_predictor - y2)/length(y2))
#
# # but the blup doesn't, because `cov(t(X2), t(X1))` doesn't provide
# # information about `cov(y2, y1)`
# blup <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = cov(t(X2), t(X1)),
# eta = object$fit$eta)
#
# drop(crossprod(blup - y2)/length(y2))
#
# # in this case, cov(y2, y1) ~= V -this information improves the blup
# blup_V <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = V,
# eta = object$fit$eta)
#
# drop(crossprod(blup_V - y2)/length(y2))
#
# # with the null model...
# object$fit$beta <- matrix(0, nrow = nrow(object$fit$beta),
# ncol = ncol(object$fit$beta))

# # ...the linear predictor is worse...
# linear_predictor_null <- predict.plmm(object$fit,
# newX = X2,
# type = 'lp',
# lambda = object$lambda.min)
# 
# drop(crossprod(linear_predictor_null - y2)/length(y2))

# # ...but the blup is even worse, because we are left only with a noisy
# # rotation of y as an estimate
# blup_null <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = cov(t(X2), t(X1)),
# eta = object$fit$eta)
# 
# drop(crossprod(blup_null - y2)/length(y2))
# 
# # again, in this case, cov(y2, y1) ~= V, which notably improves the blup
# blup_null_V <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = V,
# eta = object$fit$eta)
# 
# drop(crossprod(blup_null_V - y2)/length(y2))

# 
# 
# # Based on this, we expect to see out of sample predictions ARE improved when
# # V can be estimated from X
# n <- 210
# p <- 1000
# p1 <- 50
# SNR <- 10

# # Generate data (X related to V)
# V <- Matrix::bdiag(matrix(0.5, n/3, n/3),
# matrix(0.5, n/3, n/3),
# matrix(0.5, n/3, n/3))
# V <- as.matrix(V)
# diag(V) <- 1
# set.seed(7)
# X1 <- t(MASS::mvrnorm(n = p, mu = rep(0, n), Sigma = V))
# j <- 1:p
# s <- rep(1, times = p)
# b <- j <= p1
# b <- b * s
# covX <- stats::var(X1) * (n - 1)/n
# beta <- b * sqrt(SNR) / sqrt(drop(crossprod(b, covX) %*% b))
# Xbeta <- X1 %*% beta
# e <- MASS::mvrnorm(n=1, rep(0, n), V)
# y1 <- Xbeta + e
# 
# # fit the model
# object <- cv.plmm(X1,
#                   y1,
#                   V,
#                   eta_star = 1,
#                   type = 'blup',
#                   penalty = "lasso",
#                   alpha = 1,
#                   standardizeX = TRUE,
#                   standardizeRtX = TRUE,
#                   rotation = TRUE,
#                   returnX = TRUE,
#                   nfolds = 5)
# 
# # generate an X2 that is similar to X1
# X2 <- X1 + rnorm(nrow(X1), 0, 0.0001)
# covX <- stats::var(X2) * (n - 1)/n
# beta <- b * sqrt(SNR) / sqrt(drop(crossprod(b, covX) %*% b))
# Xbeta <- X2 %*% beta
# e <- MASS::mvrnorm(n=1, rep(0, n), V)
# y2 <- Xbeta + e
# 
# # linear predictor works well
# linear_predictor <- predict.plmm(object$fit,
# newX = X2,
# type = 'lp',
# lambda = object$lambda.min)
# drop(crossprod(linear_predictor - y2)/length(y2))

# # and the blup improves upon the linear predictor because
# # `cov(t(X2), t(X1))` ~=`cov(y2, y1)`
# blup <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = cov(t(X2), t(X1)),
# eta = object$fit$eta)
# 
# drop(crossprod(blup - y2)/length(y2))
# 
# # in this case, knowing the true(ish) cov(y2, y1) is even more beneficial than
# # using the version estimated from `cov(t(X2), t(X1))`
# blup_V <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = V,
# eta = object$fit$eta)
# 
# drop(crossprod(blup_V - y2)/length(y2))
# 
# # with the null model...
# object$fit$beta <- matrix(0, nrow = nrow(object$fit$beta),
# ncol = ncol(object$fit$beta))
# 
# # ...the linear predictor is worse...
# linear_predictor_null <- predict.plmm(object$fit,
# newX = X2,
# type = 'lp',
# lambda = object$lambda.min)

# drop(crossprod(linear_predictor_null - y2)/length(y2))
# 
# # ...but the blup is better because we still have information about
# # `cov(y2, y1)` ~= `cov(t(X2), t(X1))`
# blup_null <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = cov(t(X2), t(X1)),
# eta = object$fit$eta)
# 
# drop(crossprod(blup_null - y2)/length(y2))
# 
# # in this case, knowing the true(ish) cov(y2, y1) gives a more drastic
# # improvement than using the version estimated from `cov(t(X2), t(X1))`
# # because we don't have any information from the linear predictor
# blup_null_V <- predict.plmm(object$fit,
# newX = X2,
# type = 'blup',
# lambda = object$lambda.min,
# no_int_X = X1,
# y = y1,
# U = object$fit$U,
# S = object$fit$S,
# covariance = V,
# eta = object$fit$eta)
# 
# drop(crossprod(blup_null_V - y2)/length(y2))
# }
# 
