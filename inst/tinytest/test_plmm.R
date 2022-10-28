# Tabitha Peter 
# Fall 2022
# NOTE: These tests are in the process of being rewritten in order to be 
#   consistent with the updates I have made to the package. 

# Test 0: Case where K = identity -------------------------------------------

# set up 
K_identity <- diag(nrow = nrow(admix$X))
lambda0 <- c(1, 0.1, 0.01, 0.001)

plmm0 <- plmm(X = admix$X,
              y = admix$y,
              K = K_identity,
              lambda = lambda0,
              penalty = "lasso")

lasso0 <- glmnet(x = admix$X,
                 y = admix$y,
                 family = "gaussian",
                 lambda = lambda0)

A0 <- as.matrix(plmm0$beta_vals[2:10, ]) 
dimnames(A0) <- NULL
B0 <- as.matrix(lasso0$beta[1:9, ]) # NB: glmnet() does not return intercept values
dimnames(B0) <- NULL

# test 1 - implementation 
tinytest::expect_equal(A0, B0, tolerance = 0.01)


# Test 1 Case where K is diagonal and lambda is 0 ---------------------------

# use 'lm' here instead of glmnet

K_diagonal <- diag(x = (rnorm(n = nrow(admix$X))^2),
                   nrow = nrow(admix$X))

plmm1 <- plmm(X = admix$X,
              y = admix$y,
              K = K_diagonal,
              # even though I am interested in lambda == 0, the arg here must be
              # a sequence 
              lambda = c(0.001, 0),
              penalty = "lasso")

print(summary(plmm1, lambda = 0))

A1 <- plmm1$beta_vals[,"0.0000"]
names(A1) <- NULL

lm1 <- lm(admix$y ~ admix$X, 
          weights = 1/diag(K_diagonal))

B1 <- lm1$coefficients
names(B1) <- NULL
B1 <- ifelse(is.na(B1), 0, B1)

# test 1: implementation 
tinytest::expect_equal(A1, B1, tolerance = 0.01)

# investigate 
head(data.frame(A1, B1)) # the signs look right, but the magnitudes are off 

# Test 2: Case where K is diagonal -------------------------------------------

lambda2 <- c(1, 0.1, 0.01, 0.001)

plmm2 <- plmm(X = admix$X,
              y = admix$y,
              K = K_diagonal,
              lambda = lambda2,
              penalty = "lasso")

lasso2 <- glmnet(x = admix$X,
                 y = admix$y,
                 family = "gaussian",
                 lambda = lambda2,
                 # weights are by INVERSE variance 
                 weights = 1/diag(K_diagonal))


A2 <- as.matrix(plmm2$beta_vals[2:10, ]) 
dimnames(A2) <- NULL
B2 <- as.matrix(lasso2$beta[1:9, ]) # NB: glmnet() does not return intercept values
dimnames(B2) <- NULL

# test 2 - implementation 
tinytest::expect_equal(A2, B2, tolerance = 0.01)

# Test 3: show that monomorphic SNPs are given beta values of 0s -------------
monomorphic <- apply(admix$X[,1:15], 2, var) == 0
monomorphic_snps <- colnames(admix$X[,1:15])[monomorphic]
# SNPs 8 and 14 are monomorphic 
fit3 <- plmm(admix$X[,1:15], admix$y)

# test 3: implementation 
tinytest::expect_equivalent(summary.plmm(fit3, quiet = T)$monomorphic_snps,
                            monomorphic_snps)


# Test 4: examine the 'untransform' function ---------------------------------

# setup lambda
lambda4 <- 0.1
nlambda <- 1


# Part 1: Compute betas on the original scale (ie. sans transformation)
small_X <- admix$X[1:10, 1:10]

small_y <- admix$y[1:10]

small_K <- relatedness_mat(small_X)

small_res1 <- ncvreg::ncvreg(X = small_X,
                             y = small_y,
                             lambda = lambda4,
                             penalty = "lasso")
og_betas <- small_res1$beta


# Part 2: Compute betas after transformations
fit4 <- plmm(X = small_X, y = small_y, lambda = lambda4)

fit4$beta_vals

# QUESTION: should I be comparing fit4$beta_vals with og_betas? 

# Test 5: case where K is simulated to have no population structure ----------

K_independent <- sim_ps_x(n = nrow(admix$X),
                          p = ncol(admix$X),
                          # supposing all individuals are independent 
                          nJ = nrow(admix$X),
                          structureX = "independent"
) |> 
  relatedness_mat()


## Appendix: Code from previous version of the package -----------------------
# library(tinytest)
# library(parallel)
# devtools::load_all(".")

# nn <- 50
# pp <- 4
# p1 <- 1
# set.seed(7)
# X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
# X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
# B <- rep(c(1, 0), times = c(p1, pp - p1))
# y <- X %*% B + X0 %*% c(1) + rnorm(nn)
# V <- tcrossprod(ncvreg::std(X))/ncol(X)

# dont_run <- FALSE


### no rotation checks ------------------------------------------------------###

# std(X)
# plmm1 <- plmm(ncvreg::std(X),
#               y,
#               V = V,
#               penalty = "lasso",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# ncv1 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "lasso", lambda = plmm1$lambda)
# glm1 <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = TRUE, lambda = plmm1$lambda)

# expect_equivalent(coef(plmm1), coef(ncv1), tol = 1e-3)
# expect_equivalent(coef(plmm1), as.matrix(coef(glm1)), tol = 1e-3)
# expect_equivalent(coef(plmm1)[-1, 1], rep(0, ncol(X))) # make sure setup lambda is working correctly

###--------------------------------------------------------------------------###

# X
# plmm2 <- plmm(X,
#               y,
#               V,
#               penalty = "lasso",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = TRUE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# ncv2 <- ncvreg::ncvreg(X, y, "gaussian", penalty = "lasso", lambda = plmm2$lambda)
# glm2 <- glmnet::glmnet(X, y, "gaussian", standardize = TRUE, lambda = plmm2$lambda)

# expect_equivalent(coef(plmm2), coef(ncv2), tol = 1e-3)
# expect_equivalent(coef(plmm2), as.matrix(coef(glm2)), tol = 1e-3)
# expect_equivalent(coef(plmm2)[-1, 1], rep(0, ncol(X)))

### make sure class raw works -----------------------------------------------###
# library(snpStats)
# data(testdata)
# Autosomes <- Autosomes[1:50, 1:100]
# Autosomes cannot have missing values
# to_impute <- which(snpStats::col.summary(Autosomes)$Call.rate < 1)
# miss <- Autosomes[, to_impute]

# imputed_mean <- apply(methods::as(miss, "numeric"), 2, function(s){
#   these <- which(is.na(s))
#   s[these] <- mean(s, na.rm = TRUE)
#   s <- snpStats::mean2g(s)
#   return(s)
# })

# Autosomes2 <- Autosomes
# Autosomes2@.Data[, to_impute] <- imputed_mean
# Autosomes <- Autosomes2

# yy <- rnorm(nrow(Autosomes))
# VV <- diag(nrow(Autosomes))

# X
# plmm_raw <- plmm(Autosomes,
#               yy,
#               VV,
#               penalty = "lasso",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = TRUE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# plmm_num <- plmm(as(Autosomes, 'numeric'),
#                  yy,
#                  VV,
#                  penalty = "lasso",
#                  alpha = 1,
#                  nlambda = 5,
#                  standardizeX = TRUE,
#                  standardizeRtX = FALSE,
#                  rotation = FALSE,
#                  returnX = FALSE)

# expect_equivalent(coef(plmm_raw), coef(plmm_num), tol = 1e-2)
# expect_equivalent(plmm_raw$lambda, plmm_num$lambda, tol = 1e-5) # are they both computing lambda sequences the same way

### no rotation + unpenalized covar checks-----------------------------------###
### first var (after int) should be included because unpenalized


### these are not passing test...pass at 1e-1 (?)
### I think this has to do with the issue in glmnet not treating the intercept/unpenalized vars equivalently
### With lambda = 0 (ols) seems fine...

# std(X)
# plmm3 <- plmm(ncvreg::std(cbind(X0, X)),
#               y,
#               V = V,
#               penalty = "lasso",
#               penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))),
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# glm3 <- glmnet::glmnet(ncvreg::std(cbind(X0, X)), y, "gaussian",
#                        standardize = FALSE, lambda = plmm3$lambda,
#                        penalty.factor = rep(c(0, 1), times = c(ncol(X0), ncol(X))))

# expect_equivalent(coef(plmm3), as.matrix(coef(glm3)), tol = 1e-3)
# expect_equivalent(coef(plmm3), as.matrix(coef(glm3)), tol = 1e-1)
# expect_equivalent(coef(plmm3)[-c(1:(1 + ncol(X0))), 1], rep(0, ncol(X)))

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
# plmm5 <- plmm(ncvreg::std(X),
#               y,
#               V,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)

# if (nrow(X) > ncol(X)){
#   ols_soln <- as.numeric(solve(t(plmm5$SUX) %*% plmm5$SUX) %*% t(plmm5$SUX) %*% plmm5$SUy)
#   expect_equivalent(coef(plmm5), ols_soln, tol = 1e-3)
# }

# X
# plmm6 <- plmm(X,
#               y,
#               V,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)

# if (nrow(X) > ncol(X)){
#   ols_soln <- lm(y ~ X)
#   expect_equivalent(as.numeric(coef(plmm6)), coef(ols_soln), tol = 1e-3)
# }

### rotation + unpenalized covar checks-----------------------------------###

# std(X)
# plmm7 <- plmm(ncvreg::std(cbind(X0, X)),
#               y,
#               V,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)

# if (nrow(cbind(X0, X)) > ncol(cbind(X0, X))){
#   ols_soln <- as.numeric(solve(t(plmm7$SUX) %*% plmm7$SUX) %*% t(plmm7$SUX) %*% plmm7$SUy)
#   expect_equivalent(coef(plmm7), ols_soln, tol = 1e-3)
# }

# X
# plmm8 <- plmm(cbind(X0, X),
#               y,
#               V,
#               penalty = "lasso",
#               alpha = 1,
#               lambda = 0, # compare to ols solutions
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)

# if (nrow(cbind(X0, X)) > ncol(cbind(X0, X))){
#   ols_soln <- lm(y ~ cbind(X0, X))
#   expect_equivalent(coef(plmm8), coef(ols_soln), tol = 1e-3)
# }

### mcp checks --------------------------------------------------------------###
# must use std(X), no rotation - no other package for comparison

# plmm9 <- plmm(ncvreg::std(X),
#               y,
#               V,
#               penalty = "MCP",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# ncv9 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "MCP", lambda = plmm9$lambda)
# expect_equivalent(coef(plmm9), coef(ncv9), tol = 1e-12)


# enet
# plmm10 <- plmm(ncvreg::std(X),
#               y,
#               V,
#               penalty = "MCP",
#               alpha = 0.5,
#               nlambda = 5,
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# ncv10 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "MCP", alpha = 0.5, lambda = plmm10$lambda)
# expect_equivalent(coef(plmm10), coef(ncv10), tol = 1e-12)

### scad checks -------------------------------------------------------------###
# must use std(X), no rotation

# plmm11 <- plmm(ncvreg::std(X),
#               y,
#               V,
#               penalty = "SCAD",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = FALSE,
#               standardizeRtX = FALSE,
#               rotation = FALSE,
#               returnX = FALSE)

# ncv11 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "SCAD", lambda = plmm11$lambda)
# expect_equivalent(coef(plmm11), coef(ncv11), tol = 1e-12)

# enet
# plmm12 <- plmm(ncvreg::std(X),
#                y,
#                V,
#                penalty = "SCAD",
#                alpha = 0.5,
#                nlambda = 5,
#                standardizeX = FALSE,
#                standardizeRtX = FALSE,
#                rotation = FALSE,
#                returnX = FALSE)

# ncv12 <- ncvreg::ncvreg(ncvreg::std(X), y, "gaussian", penalty = "SCAD", alpha = 0.5, lambda = plmm12$lambda)
# expect_equivalent(coef(plmm12), coef(ncv12), tol = 1e-12)



### This doesn't work in general for glmnet - it's not the same for a manual int with int = F
# https://stackoverflow.com/questions/49495494/glmnet-is-different-with-intercept-true-compared-to-intercept-false-and-with-pen=
# if (dont_run){
#   fit <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = FALSE, intercept = TRUE, nlambda = 5)
#   fit1 <- glmnet::glmnet(cbind(1, ncvreg::std(X)), y, "gaussian", standardize = FALSE, intercept = FALSE,
#                          penalty.factor = c(0, rep(1, ncol(X))), lambda = fit$lambda)
#   coef(fit)
#   coef(fit1)
# 
#   fit <- glmnet::glmnet(ncvreg::std(X), y, "gaussian", standardize = FALSE, intercept = TRUE, lambda = 0)
#   fit1 <- glmnet::glmnet(cbind(1, ncvreg::std(X)), y, "gaussian", standardize = FALSE, intercept = FALSE,
#                          penalty.factor = c(0, rep(1, ncol(X))), lambda = 0)
#   coef(fit)
#   coef(fit1)
# }
### can't do a rotated and standardized version check - standardization occurs
### for the unrotated version of X, glmnet would give the standardization for the
### rotated version, which will not be the same.
### Using the unrotated test to check that standardization is working as it should be


