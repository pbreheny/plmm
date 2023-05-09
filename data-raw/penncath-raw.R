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
# For now, I will store this in here (data-raw/)

# Step 2: Use PLINK to arbitrarily filter out >95% of the SNPs in the data.
# This should reduce the 800K+ SNPs down to a manageable number 
# The new files will be called "penncath_lite" (the smallest version) and "_mid" (a mid-sized version)
# NOTE: I will use PLINK 1.9 for this step
system("plink --bfile penncath_og --thin 0.005 --make-bed --out penncath_lite")
system("plink --bfile penncath_og --thin 0.05 --make-bed --out penncath_mid")

# Step 3: Make sure the new files look right 
penncath_lite_fam <- fread('penncath_lite.fam')
penncath_lite_bim <- fread('penncath_lite.bim')

penncath_mid_fam <- fread('penncath_mid.fam')
penncath_mid_bim <- fread('penncath_mid.bim')

penncath_og_fam <- fread('penncath_og.fam')
penncath_og_bim <- fread('penncath_og.bim')

# Step 4: Save these files in a way that is consistent with the best practices 
# outlined in R-pkgs 2nd edition: 
# gzip the _lite and _mid versions 


# save the gzipped files in inst/extdata/


# create .rda object for penncath_mid
