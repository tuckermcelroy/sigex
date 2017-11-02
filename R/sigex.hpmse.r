sigex.hpmse <- function(N,param,mdl,trendcyclecomp,sigcomps,grid,snr)
{

	###############################
	#   sigex.hpmse
	#	by Tucker McElroy	
	#
	#	Computes signal extraction mse arising from HP filtering
	#		of trend-cycle, with HP q = snr, for trend and cycle each.
	#		Presumes that signal involves at most 2 trend differences
	#		param must be in format yielded by sigex.default
	#		trendcyclecomp is the (single) index of the trend-cycle component
	#		sigcomps provides indices of a desired component that
	#			is disjoint from trend-cycle, so that MSEs of 
	#			trend+sigcomps and cycle+sigcomps are computed
	#		 pass in sigcomps = NULL to just get trend and cycle MSEs
	#		grid is desired resolution of frequencies 
	#			(take to be a large prime, e.g., 997)
	#     Output is list of NxN matrix, for trend and cycle respectively
	#
	#################################

quad <- function(z)
{
	# quadrature of z 
	len <- length(z)
	out <- (z[1]+z[len])/2 + sum(z[2:(len-1)])
	return(out/len)
}

	lambda <- pi*seq(0,grid)/grid
	f.sig <- t(rep(0,grid+1) %x% diag(N))
	for(i in 1:length(mdl[[3]]))
	{
		L.par <- param[[1]][[i]]
		D.par <- param[[2]][[i]]
		delta <- sigex.delta(mdl,i)
		f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],N,delta,grid)
		if(i %in% trendcyclecomp) f.sig <- f.sig + matrix(f.comp,nrow=N)
	}
	frf.wksig <- sigex.wkmse(data,param,mdl,trendcyclecomp,grid)
	frf.wksig <- matrix(frf.wksig,nrow=N)

	# get first and second terms of MSE 

	allcomps <- seq(1,length(mdl[[3]]))
	noisecomps <- allcomps[!allcomps %in% trendcyclecomp]
	delta <- sigex.delta(mdl,noisecomps)
	d <- length(delta) - 1 - sum(abs(Arg(polyroot(delta))) > 10^(-8))
	integrand1 <- (1/snr)*(2 - 2*cos(lambda))^(2-d)/(1 + (1/snr)*(2 - 2*cos(lambda))^2)^2
	integrand1 <- t(integrand1 %x% matrix(1,N,N)) * f.sig
	integrand2 <- (1/snr)*(2 - 2*cos(lambda))^(2)/(1 + (1/snr)*(2 - 2*cos(lambda))^2)^2
	integrand2 <- t(integrand2 %x% matrix(1,N,N)) * frf.wksig
	integrand <- integrand1 - integrand2
	integrand <- array(integrand,c(N,N,(grid+1)))
	mse <- Re(apply(integrand,c(1,2),quad))
	
	# get third and fourth terms of MSE 
 
	seminoisecomps <- allcomps[!allcomps %in% c(trendcyclecomp,sigcomps)]
	if(length(sigcomps)==0) { frf.wksig <- matrix(0,nrow=N,ncol=(N*(grid+1))) } else {
		frf.wksig <- sigex.wkmse(data,param,mdl,sigcomps,grid)
		frf.wksig <- matrix(frf.wksig,nrow=N) }
	if(length(seminoisecomps)==0) { frf.wknoise <- matrix(0,nrow=N,ncol=N*(grid+1)) } else {
		frf.wknoise <- sigex.wkmse(data,param,mdl,seminoisecomps,grid)
		frf.wknoise <- matrix(frf.wknoise,nrow=N) }
	integrand1 <- (1 + (1/snr)*(2 - 2*cos(lambda))^2)^(-1)
	integrand2 <- (1/snr)*(2 - 2*cos(lambda))^2*(1 + (1/snr)*(2 - 2*cos(lambda))^2)^(-1)

	integrand <- t(integrand2 %x% matrix(1,N,N)) * frf.wksig +
				t(integrand1 %x% matrix(1,N,N)) * frf.wknoise
	integrand <- array(integrand,c(N,N,(grid+1)))
	mse.trend <- mse + Re(apply(integrand,c(1,2),quad))

	integrand <- t(integrand1 %x% matrix(1,N,N)) * frf.wksig +
				t(integrand2 %x% matrix(1,N,N)) * frf.wknoise
	integrand <- array(integrand,c(N,N,(grid+1)))
	mse.cycle <- mse + Re(apply(integrand,c(1,2),quad))

	return(list(mse.trend,mse.cycle))
}
