# TKP 
# March 2024
# Objective: create a PLINK version of the 'admix' data

# library(plmm) # or devtools::load_all()

# step 1: create .ped and .map files by hand 
set.seed(123)
admix_ped <- data.frame(
  FID = 1:nrow(admix$X),
  IID = 1:nrow(admix$X),
  FATHER = 0,
  MOTHER = 0,
  SEX = sample(1:2, nrow(admix$X), replace = T)
)


set.seed(123)
admix_map <- data.frame(
  CHR = sample(1:22, ncol(admix$X), replace = T),
  SNP = colnames(admix$X),
  CM = 0,
  BP = sample(15518262:218890256, ncol(admix$X)) # fake BP values
)

# step 2: create a .txt file with our outcome 


# step 3: convert .ped/.map to .bed/.bim/.fam, and make 'y' our outcome 


