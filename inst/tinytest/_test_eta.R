# test 0 ----
# run what's below in the master branch and the dev branch 
foo0 <- plmm(admix$X, admix$y) 
foo0$eta 
# compare estimated eta value between branches 

# test 1 -----
cont_phen <- get_data("~/drives/cleft/data/phs000774/qc/subsets/cont_phen/mcw_maf_lite")
X1 <- cont_phen$X[1:100, 1:1000]
y1 <- rnorm(n = nrow(X))
foo1 <- plmm(X, y1) 
foo1$eta 

# test 2 -----
hat_eta <- hat_beta0 <-rep(NA_integer_, 100)
for(i in 1:100){
  res <- test_eta_estimation(sig_s = 5, sig_eps = 3, beta0 = 0.3,
                             K = relatedness_mat(X1))
  hat_eta[i] <- res$hat_eta
  hat_beta0[i] <- res$hat_beta0
}

summary(hat_eta); 5/8
summary(hat_beta0) # very good estimate

# test 3 ---- 
# look at K when correlations are low 
hat_eta <- rep(NA_integer_, 100)
for(i in 1:100){
  K <- diag(x = (rnorm(n = 100)^2),
                     nrow = 100)
 res <- test_eta_estimation(sig_s = 5, sig_eps = 3, beta0 = 0.3,
                            K = K)
 
 hat_eta[i] <- res$hat_eta
 hat_beta0[i] <- res$hat_beta0
}

summary(hat_eta); 5/8
summary(hat_beta0)

# test 4 -----
# Need help here...
sig_s <- 5
sig_eps <- 3
true_eta <- sig_s/(sig_s + sig_eps)

K <- tcrossprod(X)*1/ncol(X)
eigen_K <- eigen(K)
# w <- (true_eta*eigen_K$values + (1 - true_eta))^(-1/2)
# wUt <- sweep(t(eigen_K$vectors), MARGIN =  1, STATS = w, FUN = "*")
# rot_X <- wUt%*%X
Sig <- true_eta*K + (1 - true_eta)*diag(nrow(K))
y <- mvtnorm::rmvnorm(n = 1,
                      mean = X%*%true_beta,
                      sigma = Sig) |> drop()
foo4 <- plmm(X, y)
foo4$eta
