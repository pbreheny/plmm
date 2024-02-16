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

# look at lippert 2011 method ---------------------------------------
# (here, y is assumed to be mean 0 on rotated scale...)


