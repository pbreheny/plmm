#' A function to generate correlated data for testing purposes
#'
#' @param n Number of observations 
#' @param p Number of features 
#' @param s Number of significant (truly associated) features
#' @param gamma Magnitude of greatest group mean confounding 
#' @param beta Magnitude of true coefficient values correpsonding to sigificant features
#' @param B Number of groups (e.g., batches, families)
#'
#' @return a list with 'toy' data 
#' @keywords internal
#'
toy_corr_dat <- function(n=100, p=256, s=4, gamma=6, beta=2, B=20) {
  mu <- matrix(rnorm(B*p), B, p)
  z <- rep(1:B, each=n/B)
  X <- matrix(rnorm(n*p), n, p) + mu[z,] |>
    ncvreg::std()
  b <- rep(c(beta, 0), c(s, p-s))
  g <- seq(-gamma, gamma, length=B)
  y <- X %*% b + g[z]
  Z <- model.matrix(~0+factor(z))
  list(y=y, X=X, beta=b, Z=Z, gamma=g, mu=mu, id=z)
}