
#' Generate multivariate trend component
#'
#' generates a component that is first difference stationary
#'
#' @param n length of series
#' @param Ndim dimensionality of each observation
#' @param Sig covariance martrix white noise component
#' @param burn burn in (defaults to 1000)
#'
#' @return Ndim x n matrix of observations
#' @export
#'

gen_trendComp = function(n, Ndim, Sig, burn=1000){
  N = n+burn
  w = t(rmvnorm(n = N, mean = rep(0,Ndim), sigma = Sig))
  s = w[,-1] - w[,-N]
  return(t(s[, (N-n):(N-1)]))
}
