# TKP 
# September 2023 

# Test 0: For diagonal K, does our LMM match LM? ----------------------------------
lmm0 <- lmm(X = admix$X, y = admix$y, K = diag(nrow(admix$X)))
lm0 <- lm(admix$y ~ (admix$X))
# NB: lm() and lmm() handle constant features differently
ns <- which(lmm0$beta_vals != 0)
tinytest::expect_equivalent(lmm0$beta_vals[ns], lm0$coefficients[ns])

# Test 1: does LMM match PLMM for same data set? -------------------------------
# use this test to see if 'untransform()' is working like it should
lmm1 <- lmm(X = pedigree$X, y = pedigree$clinical$y, K = pedigree$K)
plmm1 <- plmm(X = pedigree$X, y = pedigree$clinical$y, K = pedigree$K)
# where lambda is smallest, these values should be close:
tinytest::expect_equivalent(lmm1$beta_vals,
                            plmm1$beta_vals[,which.min(plmm1$lambda)],
                            tolerance = 0.001)

