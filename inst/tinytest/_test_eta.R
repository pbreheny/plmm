cont_phen <- get_data("~/drives/cleft/data/phs000774/qc/subsets/cont_phen/mcw_maf_lite")

# X <- admix$X # can also try this data set for testing
X <- cont_phen$X[1:100, 1:1000]

# test 1 ----
# run what's below in the master branch and the dev branch 
set.seed(1)
y <- rnorm(n = 100)
foo1 <- plmm(X, y) 
foo1$eta 
# compare estimated eta value between branches 


# test 2 -----
# Need help here -- pick up here next time 
sig_s <- 5
sig_eps <- 3
true_eta <- sig_s/(sig_s + sig_eps)

std_X <- ncvreg::std(X)
K <- tcrossprod(std_X)*(1/ncol(X))
u <- mvtnorm::rmvnorm(n = 1,
                      sigma = sig_s*K) |> drop()
eps <- mvtnorm::rmvnorm(n = 1,
                        sigma = sig_eps*diag(nrow = nrow(X))) |> drop()
true_beta <- rep(0, ncol(std_X)) # in null model, no features are included
y <- std_X%*%true_beta + u + eps

foo2 <- plmm(X, y)
foo2$eta; true_eta

# test 3 -----
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
foo3 <- plmm(X, y)
foo3$eta
