#' Generate nicely formatted lambda vec
#'
#' @param l Vector of lambda values.
#' @return A character vector of formatted lambda value names 
#' @export
#' 
#' @examples 
#' l <- c(1, 0.1, 0.01, 0.001) 
#' lamNames(l)

## from ncvreg
lamNames <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    
    d <- min(max(d,4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}
