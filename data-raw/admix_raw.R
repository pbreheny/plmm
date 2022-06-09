## code to prepare `admix_raw` dataset goes here
# load libraries
library(dplyr)
library(penalizedLMM)
# read in data for examples 
admix_raw <- read.delim("https://s3.amazonaws.com/pbreheny-data-sets/admixture.txt")
# str(admix_raw) # includes 197 obs., 100 SNPs, and racial category 


# Simulate Y, a numeric outcome representing a phenotypic trait 
# Y will depend on race as well as some SNP information, in order to illustrate
#   population structure (as is common in GWAS studies)

# determine which 2 SNPs will be significant
set.seed(522)
signif_snps <- paste0("Snp", sample(x = 1:100, size = 2))
# determine the impact of each significant SNP on the outcome
impact <- sample(x = seq(-1, 1, length.out = 100),
                 size = length(signif_snps))
# simulate some random noise 
noise <- rnorm(n = nrow(admix_raw))
# create an outcome variable 
baseline <- noise + 
  impact[1]*admix_raw[,signif_snps[1]] +
  impact[2]*admix_raw[,signif_snps[2]] 
# NB: make Africans/African Americans have higher Y, and Japanese have lower Y
#   For illustrating typical issues of GWAS studies, it is both pedagogically 
#   convenient and representative to make the impact of race have higher 
#   magnitude than the impact of the SNPs. 
admix <- admix_raw %>%
  as_tibble() %>% 
  mutate(y = case_when(
    Race %in% c("African American", "African") ~ baseline + 2,
    Race  == "Japanese" ~ baseline - 3,
    Race == "European" ~ baseline
  ))

# create a list with the data needed for analyses 
admix <- list(X = admix %>%
                dplyr::select(-c(y, Race)) %>%
                as.matrix(),
              y = admix$y,
              race = admix$Race
)

# test this data with a couple of the penalizedLMM functions

# fit a linear model with lasso
# test1 <- penalizedLMM::lasso(X = admix$X, y = admix$y, p1 = 10)
# see if the lasso model finds the truly significant SNPs
# signif_snps
# test1$coef[test1$coef != 0] 
# compare to a penalized lmm 
# RRM <- genRelatednessMat(X = scale(admix$X))
# test2 <- plmm_lasso(X = admix$X,
#                     y = admix$y,
#                     p1 = 10,
#                     V = RRM)
# test2$coef[test2$coef != 0]

# NB: for a different type of 'story' for illustration, run the simulation with 
# set.seed(123) - here, lasso finds neither of the significant SNPs, and plmm
#   finds only one (of the two). 
usethis::use_data(admix_raw, overwrite = TRUE)
