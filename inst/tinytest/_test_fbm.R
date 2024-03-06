# March 2024 
# Objective: work through development of FBM method; will use `admix` data 

# create rds object with file-backed, standardized X (one time only)
process_plink(data_dir = "data-raw", prefix = "admix") 

# examine object 
admix_rds <- readRDS(file = "data-raw/admix.rds")
str(admix_rds)

# attempt to fit filebacked model 
plmm(X = "data-raw/admix", trace = T)

# compare with in-memory model 
plmm(X = admix$X, y = admix$y, std_needed = T)
