## code to prepare 'oav' data is here
# this data comes from the study in Aline Petrin's lab
# load libraries 
library(penalizedLMM)
library(corrplot)
# background -----------------------------------------------------------
# data come from members of an extended family. This family was first studied in 
# Richieri-Costa (2011) with an interest in the patterns of 
# oculoauriculovertebral (OAV)abnormalities observed among its members. 
# The present analysis examined genotype and phenotype data from 23 individuals 
# within this extended family, with the aim of identifying the genetic mutations 
# most strongly associated with the phenotypes of macrostomia, ptosis, and ear tags. 


# read in data ------------------------------------------------------
oav <- process_plink(data_dir = "data-raw",
                     prefix = "oav",
                     rds = TRUE,
                     impute = FALSE)
fam <- read.delim(file = "data-raw/oav.fam", sep = "", header = FALSE)

# create K, the relatedness matrix -------------------------------------
# instead of estimating relatedness, I will create K using these averages:
# https://customercare.23andme.com/hc/en-us/articles/212170668-Average-Percent-DNA-Shared-Between-Relatives
K <- diag(nrow = nrow(oav$X))
dimnames(K) <- list(fam$V2, fam$V2)
# parent-child relationships (from pedigree PPT slide)
K["11", c("13", "2", "21", "22", "23", "24", "25")] <- 0.5
K[c("13", "2", "21", "22", "23", "24", "25"), "11"] <- 0.5
K[c('13', '14'), c("15", "17", "19")] <- 0.5
K[c("15", "17", "19"), c('13', '14')] <- 0.5
K[c('6', '2'), '10'] <- 0.5
K['10', c('6', '2')] <- 0.5
K['2', c('1', '4', '5')] <- 0.5
K[c('1', '4', '5'), '2'] <- 0.5
K[c('22', '29'), '30'] <- 0.5
K['30', c('22', '29')] <- 0.5
K['39', '23'] <- 0.5
K['23', '39'] <- 0.5
K[c('40', '24'),'41'] <- 0.5
K['41', c('40', '24')] <- 0.5
# full sibling relationships
K[c('15'), c('17', '19')] <- 0.5
K[c('17', '19'), '15'] <- 0.5
K['17', '19'] <- 0.5
K['19', '17'] <- 0.5
K[c('1'), c( '4', '5')] <- 0.5
K[c( '4', '5'), '1'] <- 0.5
K['4', '5'] <- 0.5
K['5', '4'] <- 0.5
# half siblings
K['10', c('1', '4', '5')] <- 0.25
K[c('1', '4', '5'), '10'] <- 0.25
# grandparent
K['11', c('15','17', '19', '10', '4', '5', '1', '30', '39', '41')] <- 0.25
# aunt/uncle 
K[c('15', '17', '19'), c('2', '21', '22', '23', '24', '25')] <- 0.25
K[c('2', '21', '22', '23', '24', '25'), c('15', '17', '19')] <- 0.25
K[c('10', '4', '5', '1'),c('13', '21', '22', '23', '24', '25')] <- 0.25
K[c('13', '21', '22', '23', '24', '25'), c('10', '4', '5', '1')] <- 0.25
K[c('30'), c('2', '21', '13', '23', '24', '25')] <- 0.25
K[c('2', '21', '13', '23', '24', '25'), '30'] <- 0.25
K['39', c('13', '2', '21', '22', '24', '25')] <- 0.25
K[c('13', '2', '21', '22', '24', '25'), '39'] <- 0.25
K['41', c('13', '2', '21', '22', '23', '25')] <- 0.25
K[c('13', '2', '21', '22', '23', '25'), '41'] <- 0.25
# 1st cousin
K[c('15', '17', '19'), c('10','4', '5', '1', '30', '39', '41')] <- 0.125
K[c('10','4', '5', '1', '30', '39', '41'), c('15', '17', '19')] <- 0.125
K[c('10','4', '5', '1'), c('30', '39', '41')] <- 0.125
K[c('30', '39', '41'), c('10','4', '5', '1')] <- 0.125
K['30', c('39', '41')] <- 0.125
K[c('39', '41'), '30'] <- 0.125
K['39', '41'] <- 0.125
K['41', '39'] <- 0.125


Kplot <- corrplot(K,
                  tl.col = "grey50",
                  is.corr = FALSE,
                  col = COL1('Purples', 200)) # looks right


# create rda object -------------------------
oav$K <- K
oav$y <- fam$V6

usethis::use_data(oav, overwrite = TRUE)
