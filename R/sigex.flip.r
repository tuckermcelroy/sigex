sigex.flip <- function(gamma,delta)
{

	###########################
	#	sigex.flip
	#		by Tucker McElroy
	#
	#   Computes Delta' * Gamma^{-1}
	#	Delta is (lag)x(lag+d) dimensional diff matrix
	#	 with rows given by coefficients of delta, a polynomial
	#	 of degree d.
	#	Gamma is the lag dimensional Toeplitz matrix from
	#	 gamma, which are autocovariances of length lag
	#	 (so acf at lags zero through lag-1).
	#	(Computation uses nested recursions.)
	#
	##############################

	d <- length(delta)-1
	Tdiff <- length(gamma)
	Delta <- matrix(0,nrow=Tdiff,ncol=(Tdiff+d))
	for(i in 1:Tdiff) Delta[i,i:(i+d)] <- rev(delta)
	val <- t(Delta) %*% solve(toeplitz(gamma))
	return(val)
}
