spec_bauer <- function(x.acf)
{
  ## Input x.acf has lags 0 through q
  
  q <- length(x.acf)-1
  thresh <- 10^{-15}
  eps <- 1
  while(eps > thresh)
  {
    x.acf <- c(x.acf,0)
    cov.mat <- toeplitz(x.acf)
    theta <- t(chol(cov.mat))[dim(cov.mat)[1],]
    gamma <- polymult(theta,rev(theta))
    eps <- sum((c(rev(x.acf),x.acf[-1]) - gamma)^2)
    print(eps)
  }
}