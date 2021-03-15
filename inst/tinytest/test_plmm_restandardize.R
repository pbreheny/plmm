
nn <- 50
pp <- 4
p1 <- 1
set.seed(7)
X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
B <- rep(c(1, 0), times = c(p1, pp - p1))
y <- X %*% B + X0 %*% c(1) + rnorm(nn)

plmm1 <- plmm(ncvreg::std(X),
              y,
              tcrossprod(ncvreg::std(X))/ncol(X),
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(X) > ncol(X)){
  xxx <- scale_varp(plmm1$SUX[,-1])
  x0 <- plmm1$SUX[,1]
  m <- attr(xxx, "center")
  s <- attr(xxx, "scale")
  b <- coef(lm(plmm1$SUy ~ 0 + x0 + xxx))
  bb <- b[-1]/s
  b0 <- b[1]
  expect_equivalent(coef(plmm1), c(b0, bb), tol = 1e-3)
}






