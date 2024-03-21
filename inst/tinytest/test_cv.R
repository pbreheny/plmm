# Test 1 - sanity check -----------------------------------------
# 'blup' and 'response' should not give the same result 
cv1 <- cv.plmm(X = admix$X, y = admix$y, seed = 321, type = 'lp')
cv2 <- cv.plmm(X = admix$X, y = admix$y, seed = 321, type='blup')
if(identical(cv1$cve, cv2$cve))stop("BLUP and response types give same result")


# Test 2 (think of what to add here) -------------------------
# cv3 <- cv.plmm(X = admix$X, y = admix$y, k = 70, trace = T)
# cv4 <- cv.plmm(X = admix$X, y = admix$y, k = 70, type = "blup", trace = T)
# identical(cv3$cve, cv4$cve)
