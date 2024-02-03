#' A function to generate K matrices for simulation
#'
#' @param n_per_group Number of observations per group (e.g., per family, per batch, ...)
#' @param n_groups Number of groups
#' @param mu Optional vector of length n_groups with the covariance value for each group. Defaults to 0.7 for all groups.
#' @param plot Logical: should a plot of K be shown? Defaults to TRUE.
#'
#' @return A K matrix
#' @keywords internal
#' 
#' @details
#' For internal use in running tests
#' 
#'
generate_K <- function(n_per_group,
                       # plot = TRUE,
                       n_groups = 10,
                       mu = 0.7){
  
  if(length(mu) == 1){
    G <- matrix(1, n_per_group, n_per_group)*mu
    all_groups <- rep(list(G), n_groups)
  } else {
    all_groups <- vector('list', n_groups)
    for(i in 1:n_groups){
      all_groups[[i]] <- matrix(1, n_per_group, n_per_group)*mu[i]
    }
  }
  
  K <- Matrix::bdiag(all_groups) |> as.matrix()
  stopifnot(isSymmetric(K))
  

  # if(plot){
  #   # corrplot::corrplot(group_struct, is.corr = F)
  #   corrplot::corrplot(K, is.corr = F)
  # }
  
  return(K)
}

