
nn <- 50
pp <- 4
p1 <- 1
set.seed(7)
X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
B <- rep(c(1, 0), times = c(p1, pp - p1))
y <- X %*% B + X0 %*% c(1) + rnorm(nn)

### no rotation checks ------------------------------------------------------###

# std(X)
plmm1 <- plmm(ncvreg::std(X),
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv1 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "lasso", lambda = plmm1$lambda)
glm1 <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = TRUE, lambda = plmm1$lambda)

expect_equivalent(coef(plmm1), coef(ncv1), tol = 1e-3)
expect_equivalent(coef(plmm1), as.matrix(coef(glm1)), tol = 1e-3)
expect_equivalent(coef(plmm1)[-1, 1], rep(0, ncol(X))) # make sure setup lambda is working correctly

# X
plmm2 <- plmm(X,
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = TRUE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv2 <- ncvreg::ncvreg(X, y, "gaussian", penalty = "lasso", lambda = plmm2$lambda)
glm2 <- glmnet::glmnet(X, y, "gaussian", standardize = TRUE, lambda = plmm2$lambda)

expect_equivalent(coef(plmm2), coef(ncv2), tol = 1e-3)
expect_equivalent(coef(plmm2), as.matrix(coef(glm2)), tol = 1e-3)
expect_equivalent(coef(plmm2)[-1, 1], rep(0, ncol(X)))


### no rotation + unpenalized covar checks-----------------------------------###
### first var (after int) should be included because unpenalized

# std(X)
plmm3 <- plmm(ncvreg::std(cbind(X0, X)),
              y,
              X_for_K = ncvreg::std(X),
              penalty = "lasso",
              penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

ncv3 <- ncvreg::ncvreg(ncvreg::std(cbind(X0, X)), y, "gaussian",
                       penalty = "lasso", lambda = plmm3$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))
glm3 <- glmnet::glmnet(ncvreg::std(cbind(X0, X)), y, "gaussian",
                       standardize = FALSE, lambda = plmm3$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))

expect_equivalent(coef(plmm3), coef(ncv3), tol = 1e-3)
expect_equivalent(coef(plmm3), as.matrix(coef(glm3)), tol = 1e-1)
expect_equivalent(coef(plmm3)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))


# X
plmm4 <- plmm(cbind(X0, X),
              y,
              X_for_K = X,
              penalty = "lasso",
              penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
              alpha = 1,
              nlambda = 5,
              standardizeX = TRUE,
              standardizeRtX = FALSE,
              rotation = FALSE,
              returnX = FALSE)

glm4 <- glmnet::glmnet(cbind(X0, X), y, "gaussian",
                       standardize = TRUE, lambda = plmm4$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))

expect_equivalent(coef(plmm4), as.matrix(coef(glm4)), tol = 1e-1)
expect_equivalent(coef(plmm4)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))

### standardization checks --------------------------------------------------###

# standardizing SUX when rotation = FALSE is the same as standardizing original X
plmm5 <- plmm(cbind(X0, X),
              y,
              X_for_K = X,
              penalty = "lasso",
              penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = FALSE,
              returnX = FALSE)

expect_equivalent(coef(plmm4), coef(plmm5), tol = 1e-12)

# and so the same lambdas should be calculated as well
expect_equivalent(plmm4$lambda, plmm4$lambda, tol = 1e-12)


### standardization + rotation checks ---------------------------------------###
# this generates X standardized, X rotated
tmp <- plmm(ncvreg::std(X),
            y,
            penalty = "lasso",
            alpha = 1,
            lambda = 0, # compare to ols solutions
            standardizeX = FALSE,
            standardizeRtX = TRUE,
            rotation = TRUE,
            returnX = TRUE)

# then it's fit in a way that should match ols
plmm6 <- plmm(tmp$X,
              tmp$y,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              intercept = FALSE,
              penalty.factor = c(0, rep(1, ncol(X))),
              rotation = FALSE,
              returnX = FALSE)

if (nrow(X) > ncol(X)){
  ols_soln <- lm(tmp$y ~ -1 + tmp$X)
  expect_equivalent(as.numeric(coef(plmm6)), coef(ols_soln), tol = 1e-3)
}
# this should also match the results for plmm5 of test_plmm - only in the case that lambda = 0


# X
# this generates X unstandardized, X rotated
tmp <- plmm(X,
            y,
            penalty = "lasso",
            alpha = 1,
            lambda = 0, # compare to ols solutions
            standardizeX = FALSE,
            standardizeRtX = TRUE,
            rotation = TRUE,
            returnX = TRUE)

# then it's fit in a way that should match ols
plmm7 <- plmm(tmp$X,
              tmp$y,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              intercept = FALSE,
              penalty.factor = c(0, rep(1, ncol(X))),
              rotation = FALSE,
              returnX = FALSE)

if (nrow(X) > ncol(X)){
  ols_soln <- lm(tmp$y ~ -1 + tmp$X)
  expect_equivalent(as.numeric(coef(plmm7)), coef(ols_soln), tol = 1e-3)
}


### another way to check this with a standardized X - do the re-standardization on the lm side
plmm8 <- plmm(ncvreg::std(X),
              y,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = TRUE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(X) > ncol(X)){
  xxx <- ncvreg::std(plmm8$X[,-1])
  x0 <- plmm8$X[,1]
  m <- attr(xxx, "center")
  s <- attr(xxx, "scale")
  b <- coef(lm(plmm8$y ~ 0 + x0 + xxx))
  bb <- b[-1]/s
  b0 <- b[1] - crossprod(m, bb)
  expect_equivalent(coef(plmm8), c(b0, bb), tol = 1e-2)
}






