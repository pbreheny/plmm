#' mfdr: Marginal false discovery rates for PLMMs
#' Based on ncvreg::mfdr()
#' 
#' @param fit A \code{plmm} object.
#' 
#' @export
#' 
#' 
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
  rownames(df) <- lam_names(fit$lambda)
  structure(df, class=c("mfdr", "data.frame"))
}
