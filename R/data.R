#' Admix: Semi-simulated SNP data 
#'
#' A dataset containing the 100 SNPs, a demographic variable representing race,
#' and a simulated outcome
#' @format A list with 3 components
#' \describe{
#'   \item{X}{SNP matrix (197 observations of 100 SNPs)}
#'   \item{y}{vector of simulated (continuous) outcomes}
#'   \item{race}{vector with racial group categorization: # 0 = African, 1 = African American, 2 = European, 3 = Japanese}
#' }
#' @source \url{https://hastie.su.domains/CASI/}
"admix"

#'  Pedigree: mock genotype data from a family-based design. 
#'
#' A mock dataset containing the genotypes, family relationships, and continuous phenotype
#' This dataset was inspired by one of T. Peter's collaborative projects involving
#' an analysis of congenital disorders in family-based data. The outcome indicated 
#' the severity of the phenotype, with higher values -> more severe. 
#' @format ## `pedigree`
#' A list with 3 components:
#'  * X: Design matrix representing 23 family members and 5 genes
#'  * K: Matrix representing the family members' relationships as expected proportions of genetic overlap. 
#'    * Values of 0.5 represent parent-child and sibling relationships; see data-raw/pedigree.R for details. 
#'    * Note: This is NOT a correlation matrix. Use `cov2cor()` or a similar function to obtain correlations.
#'  * clinical: a data frame with 23 observations and 3 variables:
#' \describe{
#'   \item{sample.id}{integer indicating sample.id (same as IID in PLINK .fam file)}
#'    \item{sex}{the biological sex of each participant (just like PLINK, 2 = female & 1 = male)}
#'    \item{y}{numeric value corresponding to phenotype severity}
#' }
#' 
"pedigree"