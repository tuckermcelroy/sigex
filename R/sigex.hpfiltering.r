sigex.hpfiltering <- function(mdl,data,trendcyclecomp,sigcomps,psi,snr,grid,window,trendFlag)
{

	###########################
	#	sigex.hpfiltering
	#		by Tucker McElroy
	#
	#	Takes data as TxN,
	#	and applies HP*WK filter with snr parameter   to
	#	each component, where the filter has been truncated to 
	#	length 2*window + 1.  The output will be HP*WK filter applied to data,
	#	  which is forecast and aftcast extended by window units, covering time points
 	#	  1-window ..., T+window.  This gives trend and cycle at times 1,...,T.
	#	Actual window used is determined by window parameter, and additional length
	#	 determined by HP snr.
	#	Then adds signal to this estimate,
	#	  yielding trend+signal and cycle+signal.  Also gives MSEs.
	#	If not signal desired, put signal = NULL and sigcomps = NULL
	#	trendFlag = TRUE for trend+signal, else get cycle+signal
	#
	####################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	param <- sigex.psi2par(psi,mdl,data)

	# determine approporiate window size based on input window and snr
	len <- window + ceiling(4*pi/acos(1 - sqrt(hp.snr)/2))
#	hp.coefs <- hpFilt(snr,window)[[1]]
	
	hpwk.out <- sigex.hpwk(data,param,mdl,trendcyclecomp,grid,len,snr)
	hpwk.filter <- hpwk.out[[1]]
	ihpwk.filter <- hpwk.out[[2]]

	hp.mse <- sigex.hpmse(N,param,mdl,trendcyclecomp,sigcomps,grid,snr)
	if(trendFlag) 
	{ 
		hp.mse <- hp.mse[[1]] 
		psi.filter <- hpwk.filter
	} else 
	{ 
		hp.mse <- hp.mse[[2]] 
		psi.filter <- ihpwk.filter
	}	

	if(length(sigcomps) > 0)
	{
		extract.signal <- sigex.wkextract(psi,mdl,data,sigcomps,grid,window,0)
	}

	if((len) > 0) { 
		leads <- c(-rev(seq(0,len-1)),seq(1,T),seq(T+1,T+len))
	} else { leads <- seq(1,T) }
	data.ext <- t(sigex.cast(psi,mdl,data,leads,FALSE))

#	n <- window
#	hp.coefs <- hpFilt(snr,n)[[1]]
#	hp.coefs <- c(rev(hp.coefs),hp.coefs[-1])

	hp.signal <- NULL
	upp <- NULL
	low <- NULL
	for(j in 1:N) 
	{
		output.j <- rep(0,T)
		for(k in 1:N)
		{
			output.k <- filter(data.ext[,k],psi.filter[j,k,],
				method="convolution",sides=2)
			output.k <- output.k[(len+1):(len+T)]
			output.j <- output.j + output.k
		}
		if(length(sigcomps) > 0) { output.j <- output.j + extract.signal[[1]][,j] }
		hp.signal <- cbind(hp.signal,output.j)
	 	mse <- hp.mse[j,j]
	 	upp <- cbind(upp,output.j + 2*sqrt(mse))
		low <- cbind(low,output.j - 2*sqrt(mse))
	}	
	
	return(list(hp.signal,upp,low))
}



