
#' Generate multivariate seasonal monthly component
#'
#' generates a component that is monthly (12) difference stationary
#'
#' @param n length of series
#' @param Phi seasonal autoregressive parameter
#' @param Ndim dimensionality of each observation
#' @param Sig covariance martrix white noise component
#' @param burn burn in (defaults to 1000)
#'
#' @return Ndim x n matrix of observations
#' @export
#'

gen_monthComp = function(n, Phi, Ndim, Sig, burn=1000){
  N = n+burn
  w = t(rmvnorm(n = N, mean = rep(0,Ndim), sigma = Sig))
  s = w[,-(1:11)] - w[,-((N-10):N)]
  return(t(s[, (N-n-11):(N-12)]))
}
