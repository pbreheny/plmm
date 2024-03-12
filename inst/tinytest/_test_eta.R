# Testing different ways to estimate eta 

# Tests using the admix data-----------
# prepare a high dim version of admix data as well
keep50 <- sample(x = 1:nrow(admix$X), size = 50)
admix_hd <- list(
  X = admix$X[keep50, ],
  y = admix$y[keep50],
  race = admix$race[keep50]
)

K <- relatedness_mat(admix$X)
K_hd <- relatedness_mat(admix_hd$X)

## compare eta estimates from simulated outcome ---------------
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

## compare using the real outcome --------
fit1 <- plmm(X = admix$X, y = admix$y); fit1$eta
fit2 <- plmm(X = admix$X, y = admix$y, lippert_eta = TRUE); fit2$eta

# high-dim case 
fit3 <- plmm(X = admix_hd$X, y = admix_hd$y); fit3$eta
fit4 <- plmm(X = admix_hd$X, y = admix_hd$y, lippert_eta = TRUE); fit4$eta


# construct an example with skewed y -------------------------------
K <- relatedness_mat(admix$X)
hat_eta <- rep(NA_integer_, 100)
y_skew <- matrix(nrow = 100, ncol = nrow(K))
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation(sig_s = 2,
                             sig_eps = 1,
                             K = K,
                             y_dist = "skewed",
                             return_y = TRUE)
  
  y_skew[i,] <- res$y
  hat_eta[i] <- res$eta
  setTxtProgressBar(pb, i)
}

# hist(y_skew[2,])
# summary(hat_eta)

# Look at a toy K matrix --------------------------------------------
# TODO: is this a sensible way to test eta estimation?? Still thinking about this...
# create a K matrix with clear structure 
# group_struct <- matrix(0.3, nrow = 3, ncol = 3)
# diag(group_struct) <- 0.9
# K <- Matrix::bdiag(group_struct, group_struct, group_struct,
#                    group_struct, group_struct, group_struct,
#                    group_struct, group_struct, group_struct,
#                    group_struct, group_struct, group_struct)
# corrplot::corrplot(as.matrix(K), is.corr = F)
# # corrplot::corrplot(as.matrix(K) |> cov2cor())
# hat_eta <- rep(NA_integer_, 100)
# pb <- txtProgressBar(0, 100, style = 3)
# for(i in 1:100){
#   res <- test_eta_estimation(sig_s = 2,
#                              sig_eps = 1,
#                              K = as.matrix(K))
#   hat_eta[i] <- res
#   setTxtProgressBar(pb, i)
# }
# 
# summary(hat_eta); boxplot(hat_eta)
# 
