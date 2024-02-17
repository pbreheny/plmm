#' a function to rotate data (or 'precondition' data) onto the transformed scale
#'
rotate <- function(s, U, eta, y = NULL, X = NULL){
  
  # we will need the sum of the nonzero values from the diagonal matrix of weights
  w2 <- ((eta*s) + (1 - eta))
  w <- w2^(-1/2)
  # create projection matrix 
  wUt <- sweep(x = t(U), MARGIN = 1, STATS = w, FUN = "*")
  
  # preconditioning
  if (!is.null(y)){
    wUt%*%y
  } else if (!is.null(X)){
    wUt%*%X
  }
}