#' multiply two power series coeficient vectors
#'
#' @param a power series coefficient vector
#' @param b power series coefficient vector
#'
#' @return cofficient vector of product
#' @export
#'

polymult <- function(a,b) {
  bb <- c(b,rep(0,length(a)-1))
  B <- toeplitz(bb)
  B[lower.tri(B)] <- 0
  aa <- rev(c(a,rep(0,length(b)-1)))
  prod <- B %*% matrix(aa,length(aa),1)
  return(rev(prod[,1]))
}

