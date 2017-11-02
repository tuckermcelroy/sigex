sigex.signalextract <- function(data,psi,mdl,sigcomps,window,grid)
{
	###############################
	#   sigex.signalextract
	#	by Tucker McElroy	
	#
	#	Combines sigex.signal and sigex.extract using a windowing 
	#		method for longer series.  
	#		Utilizes wk filter for middle of sample	
	#		psi must be in standard format 
	#		sigcomps provides indices of the desired components
	#		window controls approximation at the boundary
	#		grid gives resolution for WK filter
	#     Output is list of object with three length T vectors: the estimate,
	#		followed by upper and lower confidence intervals
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	
	param <- sigex.psi2par(psi,mdl,data)
	if(4*window >= T) {
		signal <- sigex.signal(data,param,mdl,sigcomps)
		extract <- sigex.extract(sample,signal,mdl,param)
		mids <- extract[[1]]
		upps <- extract[[2]]
		lows <- extract[[3]]
	} else {
		mids <- NULL
		upps <- NULL
		lows <- NULL
		
		# left sample boundary
		range <- seq(1,(2*window+1),1)
		sample <- as.matrix(data[range,])
		mdl.sample <- mdl
		for(i in 1:N) { mdl.sample[[4]][[i]] <- as.matrix(mdl[[4]][[i]][range,]) }
		signal <- sigex.signal(sample,param,mdl.sample,sigcomps)
		extract <- sigex.extract(sample,signal,mdl.sample,param)
		mids <- rbind(mids,as.matrix(extract[[1]][1:window,]))
		upps <- rbind(upps,as.matrix(extract[[2]][1:window,]))
		lows <- rbind(lows,as.matrix(extract[[3]][1:window,]))

		# middle sample 
		range <- seq((window+1),(T-window),1)					
		extract <- sigex.wkextract(psi,mdl,data,sigcomps,grid,2*window)
		mids <- rbind(mids,as.matrix(extract[[1]][range,]))
		upps <- rbind(upps,as.matrix(extract[[2]][range,]))
		lows <- rbind(lows,as.matrix(extract[[3]][range,]))

		# right sample boundary
		range <- seq((T-2*window),T,1)
		sample <- as.matrix(data[range,])
		mdl.sample <- mdl
		for(i in 1:N) { mdl.sample[[4]][[i]] <- as.matrix(mdl[[4]][[i]][range,]) }
		signal <- sigex.signal(sample,param,mdl.sample,sigcomps)
		extract <- sigex.extract(sample,signal,mdl.sample,param)
		mids <- rbind(mids,as.matrix(extract[[1]][(length(range)-window+1):length(range),]))
		upps <- rbind(upps,as.matrix(extract[[2]][(length(range)-window+1):length(range),]))
		lows <- rbind(lows,as.matrix(extract[[3]][(length(range)-window+1):length(range),]))
	}

	return(list(mids,upps,lows))
}

		