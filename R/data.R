#' Semi-simulated SNP data 
#'
#' A dataset containing the 100 SNPs, a demographic variable representing race,
#' and a simulated outcome
#' @format A list with 3 components
#' \describe{
#'   \item{X}{SNP matrix (197 observations of 100 SNPs)}
#'   \item{y}{vector of simulated (continuous) outcomes}
#'   \item{race}{vector with racial group categorization}
#' }
#' @source \url{https://hastie.su.domains/CASI/}
"admix"

#'  Oculo-auriculo-vertebral spectrum disorders: genotype data from an extended family
#'
#' A dataset containing the genotypes, family relationships, and number of OAV phenotypes present
#' @format A list with 5 components
#' \describe{
#'   \item{X}{Design matrix representing 23 family members and 8 genes}
#'   \item{fam}{A data frame analogous to a PLINK fam file}
#'   \item{map}{A data frame analogous to a PLINK bim file}
#'   \item{K}{Matrix representing the family members' relationships as proportions of genetic overlap. This is NOT a correlation matrix. Use `cov2cor()` or a similar function to obtain correlations.}
#'   \item{y}{Outcome - 'burden' of OAV as measured by number of phenotypes.}
#' }
#' @source Richieri-Costa, Antonio, and Lucilene Arilho Ribeiro. "Macrostomia, preauricular tags, and external ophthalmoplegia: a new autosomal dominant syndrome within the oculoauriculovertebral spectrum?." The Cleft palate-craniofacial journal 43.4 (2006): 429-434.
"oav"