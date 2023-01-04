# Tests for the cv.plmm function: 

# Test 1: ---------------------





# Appendix: the following are ideas/work from Anna R.'s dissertation-----------


# library(tinytest)
# library(parallel)
# devtools::load_all(".")



# nn <- 50
# pp <- 3
# p1 <- 1
# set.seed(7)
# X <- matrix(rnorm(nn * pp), nrow = nn, ncol = pp)
# X0 <- matrix(c(rnorm(nn * 1)), nrow = nn, ncol = 1)
# B <- rep(c(1, 0), times = c(p1, pp - p1))
# y <- X %*% B + X0 %*% c(1) + rnorm(nn)
# V <- matrix(rnorm(nn * pp), nn, nn)
# V <- crossprod(V) # make this symmetric
# diag(V) <- 1


# penalty = "lasso"
# alpha = 1
# nlambda = 5
# standardizeX = TRUE
# rotation = TRUE
# returnX = TRUE

### cv.plmm fit should match plmm fit ---------------------------------------###
### known similarity matrix
# plmm0 <- plmm(X,
#               y,
#               V,
#               eta_star = 1,
#               penalty = "lasso",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = TRUE,
#               standardizeRtX = TRUE,
#               rotation = TRUE,
#               returnX = TRUE)

# cv_plmm0 <- cv.plmm(X,
#                     y,
#                     V,
#                     eta_star = 1,
#                     type = 'response',
#                     penalty = "lasso",
#                     alpha = 1,
#                     lambda = plmm0$lambda,
#                     standardizeX = TRUE,
#                     standardizeRtX = TRUE,
#                     rotation = TRUE,
#                     returnX = TRUE,
#                     nfolds = 2)

# cv_plmm00 <- cv.plmm(X,
#                      y,
#                      V,
#                      eta_star = 1,
#                      type = 'individual',
#                      penalty = "lasso",
#                      alpha = 1,
#                      lambda = plmm0$lambda,
#                      standardizeX = TRUE,
#                      standardizeRtX = TRUE,
#                      rotation = TRUE,
#                      returnX = TRUE,
#                      nfolds = 2)

# cv_plmm1 <- cv.plmm(X,
#                     y,
#                     V,
#                     eta_star = 1,
#                     type = 'individual',
#                     penalty = "lasso",
#                     alpha = 1,
#                     nlambda = 5,
#                     standardizeX = TRUE,
#                     standardizeRtX = TRUE,
#                     rotation = TRUE,
#                     returnX = TRUE,
#                     nfolds = 2,
#                     seed = 7)

# expect_equivalent(coef(plmm0), coef(cv_plmm0$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(coef(plmm0), coef(cv_plmm00$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(plmm0$lambda, cv_plmm1$lambda, tol = 1e-5) # are they both computing lambda sequences the same way

### make sure class raw works -----------------------------------------------###
# library(snpStats)
# data(testdata)
# Autosomes <- Autosomes[1:50, 1:100]
# # Autosomes cannot have missing values
# to_impute <- which(snpStats::col.summary(Autosomes)$Call.rate < 1)
# miss <- Autosomes[, to_impute]
# 
# imputed_mean <- apply(methods::as(miss, "numeric"), 2, function(s){
#   these <- which(is.na(s))
#   s[these] <- mean(s, na.rm = TRUE)
#   s <- snpStats::mean2g(s)
#   return(s)
# })

# Autosomes2 <- Autosomes
# Autosomes2@.Data[, to_impute] <- imputed_mean
# Autosomes <- Autosomes2
# 
# yy <- rnorm(nrow(Autosomes))
# VV <- diag(nrow(Autosomes))

# cv_plmm_num <- cv.plmm(as(Autosomes, 'numeric'),
#                        yy,
#                        VV,
#                        eta_star = 1,
#                        type = 'individual',
#                        penalty = "lasso",
#                        alpha = 1,
#                        nlambda = 5,
#                        standardizeX = TRUE,
#                        standardizeRtX = TRUE,
#                        rotation = TRUE,
#                        returnX = TRUE,
#                        nfolds = 2,
#                        seed = 7)

# cv_plmm_raw <- cv.plmm(Autosomes,
#                        yy,
#                        VV,
#                        eta_star = 1,
#                        type = 'individual',
#                        penalty = "lasso",
#                        alpha = 1,
#                        nlambda = 5,
#                        standardizeX = TRUE,
#                        standardizeRtX = TRUE,
#                        rotation = TRUE,
#                        returnX = TRUE,
#                        nfolds = 2,
#                        seed = 7)
# 
# expect_equivalent(coef(cv_plmm_num$fit), coef(cv_plmm_raw$fit), tol = 1e-2)
# expect_equivalent(cv_plmm_num$lambda, cv_plmm_raw$lambda, tol = 1e-5) # are they both computing lambda sequences the same way

### cv.plmm fit should match plmm fit ---------------------------------------###
### estimating eta
# plmm0 <- plmm(X,
#               y,
#               V,
#               penalty = "lasso",
#               alpha = 1,
#               nlambda = 5,
#               standardizeX = TRUE,
#               standardizeRtX = FALSE,
#               rotation = TRUE,
#               returnX = TRUE)

# cv_plmm0 <- cv.plmm(X,
#                    y,
#                    V,
#                    type = 'response',
#                    penalty = "lasso",
#                    alpha = 1,
#                    lambda = plmm0$lambda,
#                    standardizeX = TRUE,
#                    standardizeRtX = FALSE,
#                    rotation = TRUE,
#                    returnX = TRUE,
#                    nfolds = 2)

# cv_plmm00 <- cv.plmm(X,
#                     y,
#                     V,
#                     type = 'individual',
#                     penalty = "lasso",
#                     alpha = 1,
#                     lambda = plmm0$lambda,
#                     standardizeX = TRUE,
#                     standardizeRtX = FALSE,
#                     rotation = TRUE,
#                     returnX = TRUE,
#                     nfolds = 2)

# cv_plmm1 <- cv.plmm(X,
#                     y,
#                     V,
#                     type = 'individual',
#                     penalty = "lasso",
#                     alpha = 1,
#                     nlambda = 5,
#                     standardizeX = TRUE,
#                     standardizeRtX = FALSE,
#                     rotation = TRUE,
#                     returnX = TRUE,
#                     nfolds = 2,
#                     seed = 7)

# expect_equivalent(coef(plmm0), coef(cv_plmm0$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(coef(plmm0), coef(cv_plmm00$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(plmm0$lambda, cv_plmm1$lambda, tol = 1e-5) # are they both computing lambda sequences the same way

# cl <- parallel::makeCluster(2)
# cv_plmm2 <- cv.plmm(X,
#                     y,
#                     V=V,
#                     type = 'response',
#                     penalty = "lasso",
#                     alpha = 1,
#                     nlambda = 5,
#                     standardizeX = TRUE,
#                     rotation = TRUE,
#                     returnX = TRUE,
#                     fold =cv_plmm1$fold,
#                     cluster = cl)
# parallel::stopCluster(cl)
#
# cv_plmm2 <- cv.plmm(X,
#                     y,
#                     type = 'response',
#                     penalty = "lasso",
#                     alpha = 1,
#                     nlambda = 5,
#                     standardizeX = TRUE,
#                     standardizeRtX = FALSE,
#                     rotation = TRUE,
#                     returnX = TRUE,
#                     fold =cv_plmm1$fold,
#                     seed = 7)
#
### make sure parallel method is equivalent to regular plmm
# expect_equivalent(coef(plmm0), coef(cv_plmm2$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(plmm0$lambda, cv_plmm2$lambda, tol = 1e-5) # are they both computing lambda sequences the same way
#
# ### make sure parallel method is equivalent to cv-plmm
# expect_equivalent(coef(cv_plmm1$fit), coef(cv_plmm2$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(cv_plmm1$lambda, cv_plmm2$lambda, tol = 1e-5) # are they both computing lambda sequences the same way
# expect_equivalent(cv_plmm1$cve, cv_plmm2$cve, tol = 1e-3) #  slight differences if in loop or parallel



