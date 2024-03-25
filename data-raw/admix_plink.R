# TKP 
# March 2024
# Objective: create a PLINK version of the 'admix' data

# library(plmm) # or devtools::load_all()
library(dplyr)
# step 1: create .ped file by hand --------------------------------
set.seed(123)
ped <- data.frame(
  FID = 1:nrow(admix$X),
  IID = 1:nrow(admix$X),
  FATHER = 0,
  MOTHER = 0,
  SEX = sample(1:2, nrow(admix$X), replace = T),
  y = admix$y
)

# need fake allele data - let's pretend all of these SNPs are T (minor) and C (major)
# minor alleles 
a1_alleles <-  admix$X |> 
  as.data.frame() |> 
  mutate(across(.cols = everything(),
                .fns = ~ case_when(
    .x == 0 ~ "C",
    .x == 1 ~ "T",
    .x == 2 ~ "T"
  ))) 

colnames(a1_alleles) <- paste0(colnames(admix$X), "_A1")
# check: snp8 should be constant 
# table(a1_alleles$Snp8_A1) # correct

# major alleles
a2_alleles <-  admix$X |> 
  as.data.frame() |> 
  mutate(across(.cols = everything(),
                .fns = ~ case_when(
                  .x == 0 ~ "C",
                  .x == 1 ~ "C", # NB: this is the difference from a1
                  .x == 2 ~ "T"
                ))) 

colnames(a2_alleles) <- paste0(colnames(admix$X), "_A2")
# check: snp8 should be constant 
# table(a2_alleles$Snp8_A2) # correct

# combine toy allele data 
alleles <- bind_cols(a1_alleles, a2_alleles) 

# PLINK expects that these columns are ordered like Snp1_A1, Snp1_A2, Snp2_A1, Snp2_A2, etc
snp_names <- sub("_A1", "", names(a1_alleles))
ordered_columns <- sapply(snp_names, function(snp) {
  paste0(snp, c("_A1", "_A2")) 
}) |> as.vector()
allele_data <- alleles |> dplyr::select(all_of(ordered_columns))
# check 
glimpse(allele_data)

admix_ped <- cbind(ped, alleles)
write.table(admix_ped,
            "data-raw/admix.ped",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
# check 
# read.delim("data-raw/admix.ped", header = F) |> View()

# step 2: create .map file by hand ----------------------------------------
admix_map <- data.frame(
  CHR = rep(1:20, each=5),
  SNP = colnames(admix$X),
  CM = 0,
  BP = 1:ncol(admix$X) # fake BP values
)
write.table(admix_map,
            "data-raw/admix.map",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)


# check 
# read.delim("data-raw/admix.map", header = F) |> View()

# step 3: convert .ped/.map to .bed/.bim/.fam, and make 'y' our outcome ----
system("plink --file data-raw/admix --make-bed --out inst/extdata/admix")
