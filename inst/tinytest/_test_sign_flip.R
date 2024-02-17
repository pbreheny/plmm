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
plot(foo1); plot(foo2)
foo1$eta; foo2$eta

# see what happens if data are high dim 
admix_50 <- sample(1:nrow(admix$X), 50) |> sort()
admix_hd <- list(
  X = admix$X[admix_50, ],
  y = admix$y[admix_50],
  race = admix$race[admix_50]
)

# since we get such different estimates for eta, let's try a simulation where 
# eta is known and then see which approach (keeping singular values = 0 or not)
# gives the accurate result 
K3 <- generate_K(4,5, mu = seq(0.7, 1.3, length.out = 5))
K3_trunc <- RSpectra::eigs(K3, k = 4)

pb <- txtProgressBar(0, 100, style = 3)
for(i in 1:100){
  hat_eta[i] <- test_eta_estimation(sig_s = 3,
                             sig_eps = 1,
                             K = K3)
  setTxtProgressBar(pb, i)
}


admix_hd$K <- relatedness_mat(admix_hd$X)
admix_hd$trunc_K <- choose_k(admix_hd$X)

foo3 <- plmm(X = admix_hd$X, y = admix_hd$y, K = admix_hd$trunc_K$svd_K$K_approx)
foo4 <- plmm(X = admix_hd$X, y = admix_hd$y, K = admix_hd$trunc_K$svd_K)

plot(foo3); plot(foo4) # notably different
foo3$eta; foo4$eta # eta estimates are polar opposites
dim(admix_hd$trunc_K$svd_K$K_approx) # U here is n x n, coming from eigen(K)
dim(admix_hd$trunc_K$svd_K$U) # U here is n x k -- eigenvalues = 0 are *excluded* 