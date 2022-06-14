#' Calculate index for which objective function ceases to be locally convex
#'
#' @param b Matrix of coefficient values.
#' @param X Design matrix.
#' @param penalty The penalty applied to the model. Either "MCP", "SCAD", or "lasso".
#' @param gamma The tuning parameter of the MCP/SCAD penalty. Default is 3 for MCP and 3.7 for SCAD.
#' @param l2 L2.
#' @param family Only "gaussian" currently supported.
#' @param penalty.factor A multiplicative factor for the penalty applied to each coefficient. If supplied, penalty.factor must be a numeric vector of length equal to the number of columns of X. The purpose of penalty.factor is to apply differential penalization if some coefficients are thought to be more likely than others to be in the model. In particular, penalty.factor can be 0, in which case the coefficient is always in the model without shrinkage.
#' @export
#' 
#' @examples 
#' 

# from ncvreg
convexMin <- function(b, X, penalty, gamma, l2, family = "gaussian", penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)

  if (penalty=="MCP") {
    if(is.null(gamma)){gamma <- 3} # TODO: verify that this is appropriate
    k <- 1/gamma
  } else if (penalty=="SCAD") {
    if(is.null(gamma)){gamma <- 3.7} # TODO: verify that this is appropriate
    k <- 1/(gamma-1)
  } else if (penalty=="lasso") {
    return(NULL)
  }
  if (l==0) return(NULL)

  val <- NULL
  for (i in 1:l) {
    A1 <- if (i==1) rep(1,p) else b[,i]==0
    if (i==l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[,i+1]==0
      U <- A1&A2
      L2 <- l2[i+1]
    }
    if (sum(!U)==0) next
    Xu <- X[,!U]
    p.. <- k*(penalty.factor[!U]!=0) - L2*penalty.factor[!U]
    if (family=="gaussian") {
      if (any(A1!=A2)) {
        eigen.min <- min(eigen(crossprod(Xu)/n - diag(p.., length(p..), length(p..)))$values)
      }
    } else {
      stop('penalizedLMM only supports family = `gaussian`!')
    }
    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}
