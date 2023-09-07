## code to prepare `pedigree` dataset goes here

# This data set is inspired by a real dataset that I encountered in my collaborative work.
# Unfortunately, I cannot include the real dataset due to IRB/HIPPA issues. 
# The pedigree_raw.rds object has the starting material for this internal data set. 
pedigree <- readRDS("data-raw/pedigree_raw.rds")

# simulate a mock relationship matrix: 
K <- diag(nrow = nrow(pedigree$X))
dimnames(K) <- list(rownames(pedigree$X), rownames(pedigree$X))

# parent-child relationships
K[c('13', '14'), c("15", "17", "19")] <- 0.5
K[c("15", "17", "19"), c('13', '14')] <- 0.5
K[c('6', '2'), '10'] <- 0.5
K['10', c('6', '2')] <- 0.5
K['2', c('1', '4', '5')] <- 0.5
K[c('1', '4', '5'), '2'] <- 0.5
# full sibling relationships
K[c('1'), c( '4', '5')] <- 0.5
K[c( '4', '5'), '1'] <- 0.5
K['4', '5'] <- 0.5
K['5', '4'] <- 0.5
# 1st cousin
K[c('3'), c('4')] <- 0.125
K[c('4'), c('3')] <- 0.125
K[c('21'), c('22')] <- 0.125
K[c('22'), c('21')] <- 0.125

# look at K 
library(corrplot)
corrplot(K,
         tl.col = "grey50",
         is.corr = FALSE,
         col = COL1('Purples', 200)) # looks right

pedigree$K <- K

usethis::use_data(pedigree, overwrite = TRUE)
