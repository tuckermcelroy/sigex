sigex.getwk <- function(data,param,mdl,sigcomps,plotit=TRUE,grid,len)
{

	###########################################
	#	sigex.getwk
	#		by Tucker McElroy
	#
	#	computes wk filter coefficients for desired signal and plots it
	#
	###########################################

	x <- t(data)
	N <- dim(x)[1]

	wk.out <- sigex.wk(data,param,mdl,sigcomps,grid,len)
	psi.comp <- wk.out[[1]]
	mse.comp <- wk.out[[2]]

	if(plotit && (N <= 3)) {
	par(mfrow = c(N,N))
	for(i in 1:N) 
	{
		for(j in 1:N)
		{
			plot(ts(psi.comp[i,j,],start=-len,frequency=1),
				xlab="Index",ylab="")
		}	
	} }

	return(list(psi.comp,mse.comp))
}


