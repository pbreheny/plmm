# Testing different ways to estimate eta 

# Tests using the estimate_eta() function derived per Lippert et al. (2011)-----

## compare to master branch -----------
# run what's below in the master branch and the dev branch 
foo0 <- plmm(admix$X, admix$y) 
foo0$eta 
# compare estimated eta value between branches 
 
## using admix data  ------------------------
K <- relatedness_mat(admix$X)
hat_eta <- rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation(sig_s = 2,
                             sig_eps = 1,
                             K = K)
  hat_eta[i] <- res
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta) # compare to true eta of sig_s/(sig_eps + sig_s)

## look at a toy K matrix -------------
# TODO: is this a sensible way to test eta estimation?? Still thinking about this...
# create a K matrix with clear structure 
group_struct <- matrix(0.3, nrow = 3, ncol = 3)
diag(group_struct) <- 0.9
K <- Matrix::bdiag(group_struct, group_struct, group_struct,
                   group_struct, group_struct, group_struct,
                   group_struct, group_struct, group_struct,
                   group_struct, group_struct, group_struct)
corrplot::corrplot(as.matrix(K), is.corr = F)
# corrplot::corrplot(as.matrix(K) |> cov2cor())
hat_eta <- rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation(sig_s = 2,
                             sig_eps = 1,
                             K = as.matrix(K))
  hat_eta[i] <- res
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta)

# Other ways to estimate eta ---------------------
## 2-param optim problem --------
K <- generate_K(4,5, mu = seq(0.7, 1.3, length.out = 5))
hat_eta <- hat_beta0 <-rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation_2param(sig_s = 3,
                             sig_eps = 1,
                             beta0 = 0.5, 
                             K = K,
                             init.vals = c(0.5, 0.1))
  hat_eta[i] <- res$hat_eta
  hat_beta0[i] <- res$hat_beta0
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta)
summary(hat_beta0) 

## 3- param optim problem -------
K <- generate_K(4,5, mu = seq(0.7, 1.3, length.out = 5))
hat_eta <- hat_r <- hat_beta0 <-rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation_3param(sig_s = 3,
                             sig_eps = 1,
                             beta0 = 0.5, 
                             K = K)
  hat_eta[i] <- res$hat_eta
  hat_r[i] <- res$hat_r
  hat_beta0[i] <- res$hat_beta0
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta)
summary(hat_r)
summary(hat_beta0) 

## allow the mean of the null model to be zero -------------
K <- generate_K(4,5, mu = seq(0.7, 1.3, length.out = 5))
hat_eta <- hat_r <- hat_beta0 <- rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation_null_mean_0(sig_s = 3,
                                    sig_eps = 1,
                                    beta0 = 0.5, 
                                    K = K)
  hat_eta[i] <- res$hat_eta
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta)

