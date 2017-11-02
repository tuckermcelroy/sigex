sigex.wkextract <- function(psi,mdl,data,sigcomps,grid,window,horizon)
{

	###########################
	#	sigex.wkextract
	#		by Tucker McElroy
	#
	#	Wrapper for WK signal extraction routine:
	#		does forecast and aftcast extention to window+horizon time units;
	#		applies WK filter for signal given by sigcomps,
	#		based on grid for Riemann integration,
	#		truncated to length 2*window + 1;
	#		generates signal at all time points 1-horizon, ..., T +horizon.
	#
      #     Output is list of object with three length T vectors: the estimate,
      #               followed by upper and lower confidence intervals
 	#		
	##############################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	param <- sigex.psi2par(psi,mdl,data)

	wk.out <- sigex.wk(data,param,mdl,sigcomps,FALSE,grid,window)
	wk.filter <- wk.out[[1]]
	wk.mse <- wk.out[[2]]
	if((window+horizon) > 0) { 
		leads <- c(-rev(seq(0,window+horizon-1)),seq(1,T),seq(T+1,T+window+horizon))
	} else { leads <- seq(1,T) }
	data.ext <- t(sigex.cast(psi,mdl,data,leads,FALSE))

	extract.sig <- NULL 
	upp <- NULL
	low <- NULL
	for(j in 1:N) 
	{
		output.j <- rep(0,T+2*horizon)
		for(k in 1:N)
		{
			output.k <- filter(data.ext[,k],wk.filter[j,k,],
				method="convolution",sides=2)
			output.k <- output.k[(window+1):(window+T+2*horizon)]
			output.j <- output.j + output.k
		}
		mse <- wk.mse[j,j]
		extract.sig <- cbind(extract.sig,output.j)
		upp <- cbind(upp,output.j + 2*sqrt(mse))
		low <- cbind(low,output.j - 2*sqrt(mse))
	}	

	return(list(extract.sig,upp,low))
}


