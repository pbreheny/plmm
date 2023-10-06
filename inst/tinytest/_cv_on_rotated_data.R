# TKP 
# Prediction with PLMM: follow-up 
# get latest changes
# devtools::install_github("areisett/penalizedLMM")
library(glmnet)

seeds <- 1:100
dat <- list()

gb_lowdim <- function(n=256, p=100, s=4, gamma=6, beta=2, B=16) {
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

gb_highdim <- function(n=100, p=256, s=4, gamma=6, beta=2, B=20) {
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

# case study -----------------------------------
l <- gb_lowdim()
cvg <- cv.glmnet(l$X, l$y)
cvp <- cv.plmm(l$X, l$y, penalty='lasso', trace = TRUE)
cvp_blup <- cv.plmm(l$X, l$y, penalty='lasso', type = "blup", returnBiasDetails = TRUE, trace = TRUE)

# MSPE
min(cvg$cvm)
min(cvp$avg_cve_over_folds)
min(cvp_blup$avg_cve_over_folds)


# MSE
crossprod(l$beta - coef(cvg)[-1]) |> drop()
crossprod(l$beta - coef(cvp)[-1]) |> drop()
crossprod(l$beta - coef(cvp_blup)[-1]) |> drop()


# simulation -----------------------------------
for (i in 1:length(seeds)){
  set.seed(seeds[i])
  dat[[i]] <- gb()
  names(dat)[i] <- paste0('sim',i)
}

# models: glment, plmm, and plmm with 'blup' 
cvg_models <- lapply(X = dat, FUN = function(l){cv.glmnet(l$X, l$y)})
cvp_models <- lapply(X = dat, FUN = function(l){cv.plmm(l$X, l$y, penalty='lasso')})
blup_models <- lapply(X = dat, FUN = function(l){cv.plmm(l$X, l$y, penalty='lasso',
                                                         type = 'blup',
                                                         returnBiasDetails = TRUE)})


# MSPE ----------------------------------------------------------------------
cvg_mspe <- lapply(X = cvg_models, FUN = function(cvg){min(cvg$cvm)}) |> unlist()
cvp_mspe <- lapply(X = cvp_models, FUN = function(cvp){min(cvp$cve)}) |> unlist()
blup_mspe <- lapply(X = blup_models, FUN = function(cvp){min(cvp$cve)}) |> unlist()

summary(cvg_mspe)
summary(cvp_mspe)
summary(blup_mspe)

data.frame(sim = seeds,
           glmnet = cvg_mspe,
           plmm = cvp_mspe,
           blup = blup_mspe) |> 
  tidyr::pivot_longer(cols = 2:4, 
                      names_to = "Model",
                      values_to = "MSPE") |> 
  ggplot2::ggplot(ggplot2::aes(x = Model, y = MSPE, color = Model)) + 
  ggplot2::geom_boxplot()



# MSE ----------------------------------------------------------------------
true_betas <- dat$sim1$beta # NB: all sims have same true beta
cvg_mse <- lapply(X = cvg_models,
                  FUN = function(cvg){crossprod(true_betas - coef(cvg)[-1]) |> drop()}) |> unlist()
cvp_mse <- lapply(X = cvp_models,
                  FUN = function(cvp){crossprod(true_betas - coef(cvp)[-1]) |> drop()}) |> unlist()
blup_mse <- lapply(X = blup_models,
                   FUN = function(cvp){crossprod(true_betas - coef(cvp)[-1]) |> drop()}) |> unlist()

data.frame(sim = seeds,
           glmnet = cvg_mse,
           plmm = cvp_mse,
           blup = blup_mse) |> 
  tidyr::pivot_longer(cols = 2:4, 
                      names_to = "Model",
                      values_to = "MSE") |> 
  ggplot2::ggplot(ggplot2::aes(x = Model, y = MSE, color = Model)) + 
  ggplot2::geom_boxplot() 

