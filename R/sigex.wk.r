sigex.wk <- function(data,param,mdl,sigcomps,plotit=TRUE,grid,len)
{

	###############################
	#   sigex.wk
	#	by Tucker McElroy	
	#
	#	Computes signal extraction filter coefficients
	#		param must be in format yielded by sigex.default
	#		sigcomps provides indices of the desired components
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
	frf.wk <- sigex.frf(data,param,mdl,sigcomps,grid)
	psi.wk <- Re(apply(frf.wk,c(1,2),quad))
	for(i in 1:len)
	{
		integrand <- t(exp(1i*lambda*i) %x% matrix(1,N,N)) * 
			matrix(frf.wk,nrow=N)
		integrand <- array(integrand,c(N,N,(grid+1)))
		next.coeff <- Re(apply(integrand,c(1,2),quad))
		psi.wk <- cbind(psi.wk,next.coeff)
		integrand <- t(exp(-1i*lambda*i) %x% matrix(1,N,N)) * 
			matrix(frf.wk,nrow=N)
		integrand <- array(integrand,c(N,N,(grid+1)))
		prev.coeff <- Re(apply(integrand,c(1,2),quad))
		psi.wk <- cbind(prev.coeff,psi.wk)
	}
	psi.wk <- array(psi.wk,c(N,N,(2*len+1)))

	### wk MSE
	frf.wksig <- sigex.wkmse(data,param,mdl,sigcomps,grid)
	mse.wksig <- Re(apply(frf.wksig,c(1,2),quad))

	if(plotit && (N <= 3)) {
	par(mfrow = c(N,N))
	for(i in 1:N) 
	{
		for(j in 1:N)
		{
			plot(ts(psi.wk[i,j,],start=-len,frequency=1),
				xlab="Index",ylab="")
		}	
	} }

	return(list(psi.wk,mse.wksig))
}


