#' A function to rotate filebacked data 
#'
#'@returns a list with 4 items: 
#'  * rot_X" `X` on the rotated scale, without re-standardizing (an FBM)
#'  * roty_y: `y` on the rotated scale (a numeric vector)
#'  * stdrot_X: `X` on the rotated scale, re-standardized (an FBM)
#'  * stdrot_X_scale: numeric vector of values used to standardize `rot_X`
#'  
#' @keywords internal
#' 
rotate_filebacked <- function(prep, ...){
 
  w <- (prep$eta * prep$s + (1 - prep$eta))^(-1/2)
  Ut <- bigstatsr::big_transpose(prep$U)
  wUt <- bigstatsr::big_apply(Ut,
                              a.FUN = function(X, ind, w, res){
                                sweep(x = X[,ind],
                                      MARGIN = 1,
                                      STATS = w,
                                      "*")},
                              a.combine = cbind,
                              w = w)
  
  # add column of 1s for intercept
  std_X_with_intcpt <- matrix(data = 0,
                              nrow = prep$std_X$nrow,
                              ncol = prep$std_X$ncol + 1) 
  std_X_with_intcpt[,1] <- rep(1, prep$std_X$nrow)
  std_X_with_intcpt <- std_X_with_intcpt |> bigstatsr::as_FBM()
  # fill in other columns with values of std_X
  bigstatsr::big_apply(prep$std_X,
                       a.FUN = function(X, ind, res){
                         res[,ind+1] <- X[,ind]
                       },
                       a.combine = cbind,
                       res = std_X_with_intcpt)
  
  # identical(std_X_with_intcpt[,2], prep$std_X[,1]) # check the preceding step
  
  # rotate X and y
  rot_X <- bigstatsr::FBM(nrow = nrow(wUt), ncol = std_X_with_intcpt$ncol)
  bigstatsr::big_apply(X = std_X_with_intcpt,
                       a.FUN = function(X,
                                        ind,
                                        wUt,
                                        res){
                         # TODO: revisit this to improve computational efficiency
                         for(i in 1:nrow(wUt)){
                           r <- wUt[i,,drop=FALSE]
                           v <- bigstatsr::big_cprodVec(X = X, y.row = r)
                           # browser()
                           res[i, ind] <- t(v)
                         }
                         
                       },
                       a.combine = rbind,
                       wUt = wUt,
                       res = rot_X)
  
  rot_y <- wUt%*%prep$y
  
  # re-scale rot_X
  rot_X_scale_info <- bigstatsr::big_scale()(rot_X)
  
  # stdrot_X_center <- rot_X_scale_info$center
  stdrot_X_scale <- rot_X_scale_info$scale[-1] # NB: here, take off 'scale' corresponding to intercept column
  stdrot_X <- big_std(X = rot_X,
                      # center = stdrot_X_center, # NB: do not re-center rotated data
                      scale = stdrot_X_scale,
                      ns = 1:ncol(rot_X)-1) # NB: the -1 -> we took off the intercept
  
  return(list(rot_X = rot_X,
              rot_y = rot_y,
              stdrot_X = stdrot_X,
              stdrot_X_scale = stdrot_X_scale))
}