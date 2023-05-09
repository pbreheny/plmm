## code to prepare `penncath_raw` dataset goes here
## NB: 'cad' stands for 'coronary artery disease', which was the subject of the 
##  research in the original reference paper: https://pubmed.ncbi.nlm.nih.gov/21239051/
## The name of the study was 'penncath', so that is how I will name the objects I create
## The 'og' names indicate the OriGinal files (which are quite large) - I got 
#   these original data files from the 'adv-gwas-tutorial' repo 

# load libraries
library(data.table)

# the objective here is to filter down the large external PLINK files into a size 
#   manageable for R package examples 

# Step 1: get original data 
# WORKING ON THIS 

# Step 2: Use PLINK to arbitrarily filter out >95% of the SNPs in the data.
# This should reduce the 800K+ SNPs down to a manageable number 
# The new files will be called "cad_lite" (the smallest version) and "cad" (a mid-sized version)
# NOTE: I will use PLINK 1.9 for this step
system("~/Desktop/plink --bfile inst/extdata/og_cad --thin 0.005 --make-bed --out inst/extdata/cad_lite")
system("~/Desktop/plink --bfile inst/extdata/og_cad --thin 0.05 --make-bed --out inst/extdata/cad_mid")

# Step 3: Make sure the new files look right 
cad_lite_fam <- fread('inst/extdata/cad_lite.fam')
cad_lite_bim <- fread('inst/extdata/cad_lite.bim')

cad_fam <- fread('inst/extdata/cad_mid.fam')
cad_bim <- fread('inst/extdata/cad_mid.bim')

