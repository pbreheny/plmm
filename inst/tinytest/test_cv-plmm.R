# #
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


# penalty = "lasso"
# alpha = 1
# nlambda = 5
# standardizeX = TRUE
# rotation = TRUE
# returnX = TRUE

### cv.plmm fit should match plmm fit ---------------------------------------###

plmm0 <- plmm(X,
              y,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              standardizeX = TRUE,
              standardizeRtX = FALSE,
              rotation = TRUE,
              returnX = TRUE)

cv_plmm0 <- cv.plmm(X,
                   y,
                   X,
                   type = 'response',
                   penalty = "lasso",
                   alpha = 1,
                   lambda = plmm0$lambda,
                   standardizeX = TRUE,
                   standardizeRtX = FALSE,
                   rotation = TRUE,
                   returnX = TRUE,
                   nfolds = 4)

cv_plmm00 <- cv.plmm(X,
                    y,
                    X,
                    type = 'individual',
                    penalty = "lasso",
                    alpha = 1,
                    lambda = plmm0$lambda,
                    standardizeX = TRUE,
                    standardizeRtX = FALSE,
                    rotation = TRUE,
                    returnX = TRUE,
                    nfolds = 4)

cv_plmm1 <- cv.plmm(X,
                    y,
                    penalty = "lasso",
                    alpha = 1,
                    nlambda = 5,
                    standardizeX = TRUE,
                    standardizeRtX = FALSE,
                    rotation = TRUE,
                    returnX = TRUE,
                    nfolds = 4,
                    seed = 7)

expect_equivalent(coef(plmm0), coef(cv_plmm0$fit), tol = 1e-5) # same overall fit?
expect_equivalent(coef(plmm0), coef(cv_plmm00$fit), tol = 1e-5) # same overall fit?
expect_equivalent(plmm0$lambda, cv_plmm1$lambda, tol = 1e-5) # are they both computing lambda sequences the same way
#
# cl <- parallel::makeCluster(2)
# cv_plmm2 <- cv.plmm(X,
#                     y,
#                     type = 'individual',
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
# ### make sure parallel method is equivalent to regular plmm
# expect_equivalent(coef(plmm0), coef(cv_plmm2$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(plmm0$lambda, cv_plmm2$lambda, tol = 1e-5) # are they both computing lambda sequences the same way
#
# ### make sure parallel method is equivalent to cv-plmm
# expect_equivalent(coef(cv_plmm1$fit), coef(cv_plmm2$fit), tol = 1e-5) # same overall fit?
# expect_equivalent(cv_plmm1$lambda, cv_plmm2$lambda, tol = 1e-5) # are they both computing lambda sequences the same way
# expect_equivalent(cv_plmm1$cve, cv_plmm2$cve, tol = 1e-3) #  slight differences if in loop or parallel
