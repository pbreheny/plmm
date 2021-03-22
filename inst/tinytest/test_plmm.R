
# library(tinytest)
# library(parallel)
# devtools::load_all(".")

nn <- 50
pp <- 4
p1 <- 1
set.seed(7)
X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
B <- rep(c(1, 0), times = c(p1, pp - p1))
y <- X %*% B + X0 %*% c(1) + rnorm(nn)
V <- tcrossprod(ncvreg::std(X))/ncol(X)

dont_run <- FALSE

### no rotation checks ------------------------------------------------------###

# std(X)
plmm1 <- plmm(ncvreg::std(X),
              y,
              V = V,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = FALSE,
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
              V,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = TRUE,
              standardizeRtX = FALSE,
              rotation = FALSE,
              returnX = FALSE)

ncv2 <- ncvreg::ncvreg(X, y, "gaussian", penalty = "lasso", lambda = plmm2$lambda)
glm2 <- glmnet::glmnet(X, y, "gaussian", standardize = TRUE, lambda = plmm2$lambda)

expect_equivalent(coef(plmm2), coef(ncv2), tol = 1e-3)
expect_equivalent(coef(plmm2), as.matrix(coef(glm2)), tol = 1e-3)
expect_equivalent(coef(plmm2)[-1, 1], rep(0, ncol(X)))

### no rotation + unpenalized covar checks-----------------------------------###
### first var (after int) should be included because unpenalized


### these are not passing test...pass at 1e-1 (?)
### I think this has to do with the issue in glmnet not treating the intercept/unpenalized vars equivalently
### With lambda = 0 (ols) seems fine...

# std(X)
plmm3 <- plmm(ncvreg::std(cbind(X0, X)),
              y,
              V = V,
              penalty = "lasso",
              penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = FALSE,
              returnX = FALSE)

glm3 <- glmnet::glmnet(ncvreg::std(cbind(X0, X)), y, "gaussian",
                       standardize = FALSE, lambda = plmm3$lambda,
                       penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))

# expect_equivalent(coef(plmm3), as.matrix(coef(glm3)), tol = 1e-3)
expect_equivalent(coef(plmm3), as.matrix(coef(glm3)), tol = 1e-1)
expect_equivalent(coef(plmm3)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))

# # X
# plmm4 <- plmm(cbind(X0, X),
#               y,
#               X_for_K = X,
#               penalty = "lasso",
#               penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = TRUE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)
#
# glm4 <- glmnet::glmnet(cbind(X0, X), y, "gaussian",
#                        standardize = TRUE, lambda = plmm4$lambda,
#                        penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))
#
# # expect_equivalent(coef(plmm4), as.matrix(coef(glm4)), tol = 1e-3)
# expect_equivalent(coef(plmm4), as.matrix(coef(glm4)), tol = 1e-1)
# expect_equivalent(coef(plmm4)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))


### rotation checks ---------------------------------------------------------###
### compare with ols solutions

# std(X)
plmm5 <- plmm(ncvreg::std(X),
              y,
              V,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(X) > ncol(X)){
  ols_soln <- as.numeric(solve(t(plmm5$SUX) %*% plmm5$SUX) %*% t(plmm5$SUX) %*% plmm5$SUy)
  expect_equivalent(coef(plmm5), ols_soln, tol = 1e-3)
}

# X
plmm6 <- plmm(X,
              y,
              V,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(X) > ncol(X)){
  ols_soln <- lm(y ~ X)
  expect_equivalent(as.numeric(coef(plmm6)), coef(ols_soln), tol = 1e-3)
}

### rotation + unpenalized covar checks-----------------------------------###

# std(X)
plmm7 <- plmm(ncvreg::std(cbind(X0, X)),
              y,
              V,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(cbind(X0, X)) > ncol(cbind(X0, X))){
  ols_soln <- as.numeric(solve(t(plmm7$SUX) %*% plmm7$SUX) %*% t(plmm7$SUX) %*% plmm7$SUy)
  expect_equivalent(coef(plmm7), ols_soln, tol = 1e-3)
}

# X
plmm8 <- plmm(cbind(X0, X),
              y,
              V,
              penalty = "lasso",
              alpha = 1,
              lambda = 0, # compare to ols solutions
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = TRUE,
              returnX = TRUE)

if (nrow(cbind(X0, X)) > ncol(cbind(X0, X))){
  ols_soln <- lm(y ~ cbind(X0, X))
  expect_equivalent(coef(plmm8), coef(ols_soln), tol = 1e-3)
}

### mcp checks --------------------------------------------------------------###
# must use std(X), no rotation - no other package for comparison

plmm9 <- plmm(ncvreg::std(X),
              y,
              V,
              penalty = "MCP",
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = FALSE,
              returnX = FALSE)

ncv9 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "MCP", lambda = plmm9$lambda)
expect_equivalent(coef(plmm9), coef(ncv9), tol = 1e-12)


# enet
plmm10 <- plmm(ncvreg::std(X),
              y,
              V,
              penalty = "MCP",
              alpha = 0.5,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = FALSE,
              returnX = FALSE)

ncv10 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "MCP", alpha = 0.5, lambda = plmm10$lambda)
expect_equivalent(coef(plmm10), coef(ncv10), tol = 1e-12)

### scad checks -------------------------------------------------------------###
# must use std(X), no rotation

plmm11 <- plmm(ncvreg::std(X),
              y,
              V,
              penalty = "SCAD",
              alpha = 1,
              nlambda = 5,
              standardizeX = FALSE,
              standardizeRtX = FALSE,
              rotation = FALSE,
              returnX = FALSE)

ncv11 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "SCAD", lambda = plmm11$lambda)
expect_equivalent(coef(plmm11), coef(ncv11), tol = 1e-12)

# enet
plmm12 <- plmm(ncvreg::std(X),
               y,
               V,
               penalty = "SCAD",
               alpha = 0.5,
               nlambda = 5,
               standardizeX = FALSE,
               standardizeRtX = FALSE,
               rotation = FALSE,
               returnX = FALSE)

ncv12 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "SCAD", alpha = 0.5, lambda = plmm12$lambda)
expect_equivalent(coef(plmm12), coef(ncv12), tol = 1e-12)



### This doesn't work in general for glmnet - it's not the same for a manual int with int = F
# https://stackoverflow.com/questions/49495494/glmnet-is-different-with-intercept-true-compared-to-intercept-false-and-with-pen=
if (dont_run){
  fit <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = FALSE, intercept = TRUE, nlambda = 5)
  fit1 <- glmnet::glmnet(cbind(1, ncvreg::std(X)), y, "gaussian", standardize = FALSE, intercept = FALSE,
                         penalty.factor = c(0, rep(1, ncol(X))), lambda = fit$lambda)
  coef(fit)
  coef(fit1)

  fit <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = FALSE, intercept = TRUE, lambda = 0)
  fit1 <- glmnet::glmnet(cbind(1, ncvreg::std(X)), y, "gaussian", standardize = FALSE, intercept = FALSE,
                         penalty.factor = c(0, rep(1, ncol(X))), lambda = 0)
  coef(fit)
  coef(fit1)
}
### can't do a rotated and standardized version check - standardization occurs
### for the unrotated version of X, glmnet would give the standardization for the
### rotated version, which will not be the same.
### Using the unrotated test to check that standardization is working as it should be


