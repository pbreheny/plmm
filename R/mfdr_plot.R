#' Plot method for mFDR 
#'
#' @param mfdr A data frame as returned by `mfdr()`
#' @param logscale Logical: should lambda values be plotted on the log scale? Defaults to TRUE.
#' @param ... Other arguments passed to `plot()`
#'
#' @export
#'
#' @examples
#' cv <- cv.plmm(admix$X, admix$y)
#' fit <- cv$fit
#' mfdr_plot(mfdr(fit), type = "l")
mfdr_plot <- function(mfdr, logscale = TRUE, ...){
  dat <- cbind(lambda = as.numeric(rownames(mfdr)),
               mfdr)
  idx <- is.finite(dat$mFDR) # look only where mFDR is finite
  if(logscale){
    ll <- log(dat$lambda)
    plot(x = ll,
         y = dat$mFDR[idx],
         xlab = expression(log(lambda)),
         ylab = "mFDR",
         xlim = c(max(ll), min(ll)),
         ...)
  } else {
    plot(x = dat$lambda,
         y = dat$mFDR,
         xlab = expression(lambda),
         ylab = "mFDR",
         ...)
  }
  
  
}