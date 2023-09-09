#' mfdr: Marginal false discovery rates for PLMMs
#' Based on [ncvreg::mfdr()]
#' 
#' @param fit A \code{plmm} object.
#' 
#' @return A data frame with a row for every value of lambda and 3 columns:
#'  * EF: The number of variables selected at the given value of lambda, averaged over the permutation fits.
#'  * S: The actual number of selected variables for the non-permuted data. 
#'  * mFDR: The estimated marginal false discovery rate (EF/S). 
#' 
#' @details
#' The function estimates the marginal false discovery rate (mFDR) for a penalized regression model. The estimate tends to be accurate in most settings, but will be slightly conservative if predictors are highly correlated
#' 
#' 
#' @export
#' 
#' @examples
#' fit <- plmm(admix$X, admix$y)
#' mfdr(fit) |> head()
mfdr <- function(fit) {
  # Setup
  if (!inherits(fit, 'plmm')) stop('"fit" must be a plmm object', call.=FALSE)
  S0 <- sum(fit$penalty.factor==0)
  S <- predict(fit, type="nvars") - S0
  
  # Call C functions
  EF <- .Call('mfdr_gaussian', fit)
  
  # Calculate rate, return
  EF <- pmin(EF - S0, S)
  mFDR <- EF/S
  mFDR[S==0] <- 0
  df <- data.frame(EF=EF, S=S, mFDR=mFDR)
  rownames(df) <- lamNames(fit$lambda)
  structure(df, class=c("mfdr", "data.frame"))
}
