# Tests for the cv.plmm function: 

# Test 1: UNDER CONSTRUCTION --------------------- 
# cv.plmm and plmm should give the same results when: 
#   (1) K is known and 
#   (2) eta = 1 and 
#   (3) penalty is lasso 

cv1 <- cv.plmm(X = admix$X,
               y = admix$y,
               K = relatedness_mat(admix$X),
               eta_star = 1,
               penalty = "lasso")

summary(cv1)

plmm1 <- plmm(X = admix$X,
              y = admix$y,
              K = relatedness_mat(admix$X),
              eta_star = 1,
              penalty = "lasso",
              alpha = 1,
              nlambda = 5,
              returnX = TRUE)

summary(plmm1, lambda = cv1$lambda.min)
plot(plmm1)
