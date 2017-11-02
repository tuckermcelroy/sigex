sigex.getfrf <- function(data,param,mdl,sigcomps,plotit=TRUE,grid)
{

	###########################################
	#	sigex.getfrf
	#		by Tucker McElroy
	#
	#	computes frf for desired signal and plots it
	#
	###########################################

	x <- t(data)
	N <- dim(x)[1]

	frf.comp <- sigex.frf(data,param,mdl,sigcomps,grid)

	if(plotit && (N <= 3)) {
	par(mfrow = c(N,N))
	for(i in 1:N) 
	{
		for(j in 1:N)
		{
			plot(ts(Re(frf.comp[i,j,]),start=0,frequency=grid),
				xlab="Cycles",ylab="")
		}
	} }

	return(frf.comp)
}


