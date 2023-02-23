#' Predict method for plmm class
#'
#' @param object An object of class \code{plmm}.
#' @param newX Design matrix used for computing predicted values if requested.
#' @param type A character argument indicating what type of prediction should be returned. Options are "response," "coefficients," "vars," "nvars," and "blup." See details. 
#' @param lambda A numeric vector of regularization parameter \code{lambda} values at which predictions are requested.
#' @param idx Vector of indices of the penalty parameter \code{lambda} at which predictions are required. By default, all indices are returned.
#' @param X Optional argument. Original design matrix (not including intercept column) from object. Required if \code{type == 'blup'} and object is too large to be returned in `plmm` object.
#' @param y Optional argument. Original continuous outcome vector from object. Required if \code{type == 'blup'}.
#' @param U Optional argument. Eigenvectors from the similarity matrix from object. Required if \code{type == 'blup'}.
#' @param S Optional argument. Eigenvalues from the similarity matrix from object. Required if \code{type == 'blup'}.
#' @param eta Optional argument. Estimated $eta$ value from object. Required if \code{type == 'blup'}.
#' @param covariance Optional argument. $q times n$ covariance matrix between new and old observations. Required if \code{type == 'blup'}.
#' @param ... Additional optional arguments
#' 
#' @rdname predict.plmm
#' @export
#'
#' @examples
#' # fit a model 
#' fit <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X))
#' 
#' # simulate new data 
#' newX <- sim_ps_x(n = nrow(admix$X), nJ = 4, p = ncol(admix$X),
#'  structureX = "independent", inbr = "heterogeneous", standardizeX = FALSE)
#'  
#'  # make predictions for all lambda values 
#'  pred1 <- predict(object = fit, newX = newX, type = "response")
#'  
#'  # make predictions for a select number of lambda values 
#'  pred2 <- predict(object = fit, newX = newX, type = "response", idx=98)
#'  
#'  # make prediction using blup 
#'  pred3 <- predict(object = fit, newX = newX, type = "blup", idx=98)
#'
#'  # compare y predictions 
#'   compare_y <- data.frame(y = admix$y, yhat_response = pred2, yhat_blup = pred3)
#'  
#'  \dontrun{
#'  fit_noX <- plmm(X = admix$X, y = admix$y, K = relatedness_mat(admix$X),
#'   returnX = FALSE)
#'  
#'  pred4 <- predict(object = fit_noX, newX = newX, type = "blup", idx=98)
#'  
#'  }
#'  
#'  
#' # make predictions when X is big
#' cad_mid <- process_plink(prefix = "cad_mid", dataDir = plink_example(path="cad_mid.fam", parent=T))
#' cad_clinical <- read.csv(plink_example(path="cad_clinical.csv"))
#' # for the sake of illustration, I use a simple mean imputation for the outcome 
#' cad_clinical$hdl_impute <- ifelse(is.na(cad_clinical$hdl), mean(cad_clinical$hdl, na.rm = T), cad_clinical$hdl)
#' fit_cad <- plmm(X = cad_mid$genotypes, y = cad_clinical$hdl_impute, k = 5)
#' cad_X <- cad_mid$genotypes
#' cad_y <- cad_clinical$hdl_impute
#' newX_cad <- sim_ps_x(n = nrow(cad_X), nJ = 4, p = ncol(cad_X),
#'  structureX = "independent", inbr = "heterogeneous", standardizeX = FALSE)
#' pred_cad <- predict(object = fit_cad, newX = newX_cad, type='blup', idx = 95, X = cad_X, y = cad_y)
#' head(data.frame(cad_y, pred_cad))
#' 
#'  
#' 

predict.plmm <- function(object, newX, type=c("response", "coefficients", "vars", "nvars", "blup"),
                           lambda, idx=1:length(object$lambda), X, y, U, S, eta, covariance, ...) {
  type <- match.arg(type)
  beta_vals <- coef.plmm(object, lambda=lambda, which=idx, drop=FALSE) # includes intercept 
  
  # addressing each type: 
  
  if (type=="coefficients") return(beta_vals)

  if (type=="nvars") return(apply(beta_vals[-1, , drop=FALSE]!=0, 2, sum)) # don't count intercept
  
  if (type=="vars") return(drop(apply(beta_vals[-1, , drop=FALSE]!=0, 2, FUN=which))) # don't count intercept
  
  Xbeta <- cbind(1, newX) %*% beta_vals
  
  if (type=="response") return(drop(Xbeta))
  
  if (type == "blup"){
    warning("The BLUP option is under development. Rely on these estimates at your own risk.")
    # case 1: the object contains all items needed for blup prediction 
    if ("X" %in% names(object) & ("SUX" %in% names(object))){
      
      # calculate covariance between new and old observations 
      covariance <- cov(t(newX), t(object$X))

      ranef <- covariance %*% object$U %*% diag((1 + object$eta * (object$S - 1))^(-1)) %*% t(object$U) %*% (object$y - cbind(1, object$X) %*% beta_vals)
      # NB: can't just use the rotated y and x here - need to scale by inverse of K, not sqrt(K)
      
      # print(eta) 
      # TODO: need to create the Xbeta object 
      blup <- Xbeta + ranef
    } 
    
    # case 2: some calculations must be done before the blup prediction
   if (!("X" %in% names(object) & ("SUX" %in% names(object)))) {
      if(missing(X) | missing(y)) stop("The design matrix is required for BLUP calculation, but is not available in the plmm object.\n This is ususally because X is large.\n Please supply the no-intercept design matrix to the X argument, and the vector of outcomes to the y argument.")
      # calculate K 
      K <- relatedness_mat(X)
      # calculate S and U 
      c(S, U) %<-% svd(K)[1:2]
      # TODO: consider how to incorporate RSpectra here 
      
      # calculate covariance between new and old observations 
      covariance <- cov(t(newX), t(X))
      
      ranef <- covariance %*% U %*% diag((1 + object$eta * (S - 1))^(-1)) %*% t(U) %*% (y - cbind(1, X) %*% beta_vals)
      # print(eta) 
      
      blup <- Xbeta + ranef
      
    }
    
    return(blup)
  }

}



# Appendix: Ideas from Anna's work 


# fit <- plmm(X = admix$X, y = admix$y, lambda = c(0.05, 0.01))
# head(predict(object = fit, newX = admix$X, type = 'response', lambda = 0.05))
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
# type = 'response',
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
# type = 'response',
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
# type = 'response',
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
# type = 'response',
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
