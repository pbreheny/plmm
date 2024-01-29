# TKP 
# Jan. 2024
# Objective: explore the 'sign flip' issue in SVD -- see if this is what
#   makes PLMM + BLUP so sensitive to truncating SVD. See Bro2007

# Idea: is base::svd() doing something to adjust for sign flipping that 
# RSpectra::svds() does not do? 
library(dplyr)
# first, a toy example ------------------------------
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

# see if 'flipping' fixes the issue
svd_K1_trunc_flipped <- flip_signs(X = K1,
                                   U = svd_K1_trunc$u,
                                   V = svd_K1_trunc$u,
                                   d = svd_K1_trunc$d)

svd_K1_trunc_flipped$U[1:5, 1:5] |> round(3)

# now, a more realistic example -------------------------------
## create an example data set ---------------------------------



## full SVD --------


## trunc SVD -------


## trunc SVD + sign check ----
