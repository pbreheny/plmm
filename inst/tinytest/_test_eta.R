# Testing different ways to estimate eta 

# Tests using the admix data-----------
 
## compare eta estimates from simulated outcome ---------------
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


# look at Lippert 2011 method (here, y is assumed to be mean 0 on rotated scale...)
lippert_hat_eta <- rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  lippert_hat_eta[i] <- lippert_test_eta_estimation(sig_s = 2,
                             sig_eps = 1,
                             K = K)
  setTxtProgressBar(pb, i)
}

summary(lippert_hat_eta); boxplot(lippert_hat_eta) # compare to true eta of sig_s/(sig_eps + sig_s)

## using the real outcome --------


# Look at a toy K matrix -------------
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

