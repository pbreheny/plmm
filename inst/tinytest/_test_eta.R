cont_phen <- get_data("~/drives/cleft/data/phs000774/qc/subsets/cont_phen/mcw_maf_lite")
X <- cont_phen$X[1:100, 1:1000]

# test 1 ----
set.seed(1)
y <- rnorm(n = 100)
foo1 <- plmm(X, y) 
foo1$eta # compare estimated eta value between branches 


# test 2 -----
# Need help here -- pick up here next time 
sig_s <- 5
sig_eps <- 3
true_eta <- sig_s/(sig_s + sig_eps)
K <- diag(100)
# Sig <- true_eta*K + (1 - true_eta)*diag(nrow(K)) # pick up here...

# test 3 -----
# Need help here...
sig_s <- 5
sig_eps <- 3
true_eta <- sig_s/(sig_s + sig_eps)
# true_beta <- sample(0:3, size = ncol(X), prob = c(0.8, 0.01, 0.05, 0.05), replace = T)
K <- tcrossprod(X)*1/ncol(X)
eigen_K <- eigen(K)
# w <- (true_eta*eigen_K$values + (1 - true_eta))^(-1/2)
# wUt <- sweep(t(eigen_K$vectors), MARGIN =  1, STATS = w, FUN = "*")
# rot_X <- wUt%*%X
Sig <- true_eta*K + (1 - true_eta)*diag(nrow(K))
y <- mvtnorm::rmvnorm(n = 1,
                      mean = X%*%true_beta,
                      sigma = Sig) |> drop()
foo3 <- plmm(X, y)
foo3$eta
