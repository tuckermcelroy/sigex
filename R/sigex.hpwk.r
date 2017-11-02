sigex.hpwk <- function(data,param,mdl,trendcyclecomp,grid,len,snr)
{

	###############################
	#   sigex.hpwk
	#	by Tucker McElroy	
	#
	#	Computes signal extraction filter coefficients
	#		for trend and cycle, by combining WK filter for trend-cycle
	#		(specified by trendcyclecomp) with HP filter of snr.
	#		param must be in format yielded by sigex.default
	#		grid is desired resolution of frequencies 
	#			(take to be a large prime, e.g., 997)
	#     Output is in array form of dimension c(N,N,2*len+1), 
	#		of real number entries, where the third dimension
	#		indexes coefficients from -len,...,0,...,len 
	#	CAUTION: take grid >> len
	#
	#################################

quad <- function(z)
{
	# quadrature of z 
	len <- length(z)
	out <- (z[1]+z[len])/2 + sum(z[2:(len-1)])
	return(out/len)
}

	x <- t(data)
	N <- dim(x)[1]
	lambda <- pi*seq(0,grid)/grid
	
	### wk filter 
	frf.wk <- sigex.frf(data,param,mdl,trendcyclecomp,grid)
	integrand.hp <- (1 + (1/snr)*(2 - 2*cos(lambda))^2)^(-1)
	integrand.ihp <- (1/snr)*(2 - 2*cos(lambda))^2*(1 + (1/snr)*(2 - 2*cos(lambda))^2)^(-1)
	frf.hpwk <- t(integrand.hp %x% matrix(1,N,N)) * matrix(frf.wk,nrow=N)
	psi.hpwk <- Re(apply(array(frf.hpwk,c(N,N,(grid+1))),c(1,2),quad))
	frf.ihpwk <- t(integrand.ihp %x% matrix(1,N,N)) * matrix(frf.wk,nrow=N)
	psi.ihpwk <- Re(apply(array(frf.ihpwk,c(N,N,(grid+1))),c(1,2),quad))
	for(i in 1:len)
	{
		integrand.hp <- t(exp(1i*lambda*i) %x% matrix(1,N,N)) * frf.hpwk
		integrand.hp <- array(integrand.hp,c(N,N,(grid+1)))
		next.coeff <- Re(apply(integrand.hp,c(1,2),quad))
		psi.hpwk <- cbind(psi.hpwk,next.coeff)
		integrand.ihp <- t(exp(1i*lambda*i) %x% matrix(1,N,N)) * frf.ihpwk
		integrand.ihp <- array(integrand.ihp,c(N,N,(grid+1)))
		next.coeff <- Re(apply(integrand.ihp,c(1,2),quad))
		psi.ihpwk <- cbind(psi.ihpwk,next.coeff)
		integrand.hp <- t(exp(-1i*lambda*i) %x% matrix(1,N,N)) * frf.hpwk
		integrand.hp <- array(integrand.hp,c(N,N,(grid+1)))
		prev.coeff <- Re(apply(integrand.hp,c(1,2),quad))
		psi.hpwk <- cbind(prev.coeff,psi.hpwk)
		integrand.ihp <- t(exp(-1i*lambda*i) %x% matrix(1,N,N)) * frf.ihpwk
		integrand.ihp <- array(integrand.ihp,c(N,N,(grid+1)))
		prev.coeff <- Re(apply(integrand.ihp,c(1,2),quad))
		psi.ihpwk <- cbind(prev.coeff,psi.ihpwk)
	}
	psi.hpwk <- array(psi.hpwk,c(N,N,(2*len+1)))
	psi.ihpwk <- array(psi.ihpwk,c(N,N,(2*len+1)))

	###  HP wk MSE
#	hp.mse <- sigex.hpmse(N,param,mdl,trendcyclecomp,NULL,grid,snr)
#	mse.hpwk <- hp.mse[[1]]
#	mse.ihpwk <- hp.mse[[2]] 	

	return(list(psi.hpwk,psi.ihpwk))
#	return(list(psi.hpwk,mse.hpwk,psi.ihpwk,mse.ihpwk))
}


