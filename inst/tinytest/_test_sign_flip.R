set.seed(123)
K1 <- generate_K(n_groups = 5, n_per_group = 4, mu = runif(5, 0.5, 2))

# compare decompositions
eigen_K1 <- eigen(K1)
U1 <- eigen_K1$vectors
s1 <- eigen_K1$values
svd_K1_trunc <- RSpectra::svds(A = K1, k = 10, nv = 10)
cat("\nU[1:5, 1:5] where U is from full decomposition")
U1[1:5, 1:5] |> round(5)
cat("\nU[1:5, 1:5] where U is from truncated decomposition")
svd_K1_trunc$u[1:5, 1:5] |> round(5)

# trying to understand how to flip signs... 
# if U has signs 'flipped', does this 'mess up' USU^T? 
admix$K <- relatedness_mat(admix$X)
admix$svd_K <- svd(admix$K)
admix$trunc_K <- choose_k(admix$X)

foo1 <- plmm(X = admix$X, y = admix$y, K = admix$trunc_K$svd_K$K_approx)
foo2 <- plmm(X = admix$X, y = admix$y, K = admix$trunc_K$svd_K)


# since we get such different estimates for eta, let's try a simulation where 
# eta is known and then see which approach (keeping singular values = 0 or not)
# gives the accurate result 
K3 <- generate_K(4,5, mu = seq(0.7, 1.3, length.out = 5))
K3_trunc <- RSpectra::eigs(K3, k = 4)
hat_eta <- hat_beta0 <-rep(NA_integer_, 100)
pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  res <- test_eta_estimation(sig_s = 3,
                             sig_eps = 1,
                             beta0 = 2,
                             K = K3)

  hat_eta[i] <- res$hat_eta
  hat_beta0[i] <- res$hat_beta0
  setTxtProgressBar(pb, i)
}

summary(hat_eta); boxplot(hat_eta)
summary(hat_beta0)

