#' a function to check for flipped signs of an SVD
#' @param X matrix whose SVD is UDV^t. SVD may be truncated. 
#' @param U matrix of left singular vectors 
#' @param d vector of singular values 
#' @param V matrix of right singular vectors 
#' @param svd_list optional list with components u, d, and v
#' 
#' @details
#' Algorithm here is inspired by Bro2007: https://www.osti.gov/servlets/purl/920802
#' 
#' @keywords internal
#' 
#' @examples
#' \dontrun{
#' svd_X <- svd(admix$X)
#' svd_X$u[1:3, 1:3]
#' flip_signs(X = admix$X, svd_list = svd_X) -> foo_mat
#' foo_mat$U[1:3, 1:3]
#' 
#' K <- relatedness_mat(admix$X)
#' K[1:4, 1:4]
#' full_svd <- svd(K)
#' full_svd$u[1:4, 1:4]
#' trunc_svd <- RSpectra::svds(K, k = 10)
#' trunc_svd$u[1:4, 1:4] # notice: signs flip (which is a 'bad sign')
#' 
#' # let's fix this: 
#' corrected_trunc_svd <- flip_signs(X = K, svd_list = trunc_svd)
#' corrected_trunc_svd$U[1:4, 1:4] 
#' }

flip_signs <- function(X, U, d, V, svd_list = NULL){
  if(!is.null(svd_list)){
    U <- svd_list$u
    d <- svd_list$d
    V <- svd_list$v
  }
  
  # initialize objects 
  terms_for_sums <- vector("list", length = length(d))
  left_signs <- rep(NA, ncol(U))
  right_signs <- rep(NA, ncol(V))
  
  # create a list of matrices, each representing UDV' for a single eigenvalue 
  Sum <- 0
  for (k in 1:length(d)){
  terms_for_sums[[k]] <- d[k]*as.matrix(tcrossprod(U[,k], V[,k]))
  Sum <- Sum + terms_for_sums[[k]]
  }

  for (k in 1:length(d)){
    F_ <- X - (Sum - terms_for_sums[[k]]) # note: "F" stands for "flip"
    
    # determine signs for U
    cp1 <- crossprod(U[,k], F_[,k]) |> drop() # CP = cross product 
    left_signs[k] <- sum(sign(cp1)*(cp1^2))
    
    # determine signs for V 
    cp2 <- crossprod(V[,k], F_[k,]) |> drop()
    right_signs[k] <- sum(sign(cp2)*(cp2^2))
    
    # determine corrected signs
    if (left_signs[k]*right_signs[k] < 0){
      if (left_signs[k] < right_signs[k]){
        left_signs[k] <- -1*left_signs[k]
      } else {
        right_signs[k] <- -1*right_signs[k]
      }
    }
    
    
  } # close for loop
  
  U <- sweep(U, 2, sign(left_signs), "*")
  V <- sweep(V, 2, sign(right_signs), "*")
  
  return(list(U = U, V = V))
  
  
}