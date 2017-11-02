sigex.precision <- function(data,param,mdl,sigcomps)
{

	###############################
	#   sigex.precision
	#	by Tucker McElroy	
	#
	#	Compares signal extraction MSE arising from
	#	 multivariate fit to implied univariate fit,
	#	 as a ratio.
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	out <- sigex.mvar2uvar(x,param,mdl)
	mdl.uni <- out[[1]]
	par.uni <- out[[2]]

	signal <- sigex.signal(x,param,mdl,sigcomps)
	signal.uni <- sigex.signal(x,par.uni,mdl.uni,sigcomps)

	mse <- matrix(diag(signal[[2]]),T,N)
	mse.uni <- matrix(diag(signal.uni[[2]]),T,N)

	return(mse/mse.uni)
}
