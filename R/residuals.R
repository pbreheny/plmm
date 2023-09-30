#' Extract residuals from a PLMM fit
#' 
#' Currently, only deviance residuals are supported.
#' 
#' @param object   Object of class `plmm`.
#' @param lambda   Values of the regularization parameter at which residuals are requested (numeric vector). For values of lambda not in the sequence of fitted models, linear interpolation is used.
#' @param which    Index of the penalty parameter at which residuals are requested (default = all indices). If `lambda` is specified, this take precedence over `which`.
#' @param drop     By default, if a single value of lambda is supplied, a vector of residuals is returned (logical; default=`TRUE`). Set `drop=FALSE` if you wish to have the function always return a matrix (see [drop()]).
#' @param unrotate Logical: should residuals be 'unrotated', i.e. transformed back to the original scale? 
#' @param ...      Not used.
#' 
#' @rdname residuals.plmm
#' 
#' @examples
#' \dontrun{
#' fit <- plmm(admix$X, admix$y)
#' residuals.plmm(fit)[1:5, 1:5]
#' head(residuals.plmm(fit, which = 50))
#' }
#' @export

residuals.plmm <- function(object,
                           lambda,
                           which=1:length(object$lambda),
                           drop=TRUE, 
                           unrotate = FALSE,
                           ...) {
  
  # Calculate matrix of residuals
  R <- matrix(nrow = nrow(object$linear.predictors), ncol = ncol(object$linear.predictors))
  for(j in 1:ncol(R)){
    R[,j] <- object$rot_y - object$linear.predictors[j]
  }
  if(unrotate){
    stop("\nThis option is still underdevelopment. Don't try to unrotate residuals yet.")
    w <- (object$eta * object$s + (1 - object$eta))^(-1/2)
    W_inv <- sweep(object$U, 2, 1/(w), "*")
    WR <- W_inv%*%R
    # TODO: if this were worth pursuing, I need to finish working this thru

  } 
  
  # Interpolate and return
  if (!missing(lambda)) {
    ind <- approx(object$lambda, seq(object$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    q <- ind %% 1
    out <- (1-q)*R[, l, drop=FALSE] + q*R[, r, drop=FALSE]
    colnames(out) <- round(lambda, 4)
  } else {
    out <- R[, which, drop=FALSE]
  }
  if (drop) return(drop(out)) else return(out)
}

# dr <- function(f, y, m) {
#   sqrt(pmax(f(y, m, rep(1, length(y))), 0)) * ((y > m) * 2 - 1)
# }