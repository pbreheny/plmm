# TKP 
# March 2024
# Objective: create a miniature version of the penncath lite data for testing purposes

# library(plmm) # I was using devtools::load_all() here

# note: I unzipped the penncath_lite data files into a temporary folder on my local machine
system("plink --bfile ~/Desktop/temp_files/penncath_lite --thin-count 500 --thin-indiv-count 100 --make-bed --out inst/extdata/penncath_mini_test")
