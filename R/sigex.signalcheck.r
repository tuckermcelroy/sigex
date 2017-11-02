sigex.signalcheck <- function(signal,param,mdl,sigcomps,lagall)
{

	###############################
	#   sigex.signalcheck
	#	by Tucker McElroy	
	#
	#	Computes a signal diagnostic, based on comparing
	#		sample autocovariances of extracted signal
	#		to the theoretical autocovariances of the 
	#		model-based signal
	#
	#################################

	x <- t(signal)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# compute differencing polynomials and matrices
	allcomps <- seq(1,length(mdl[[3]]))
	noisecomps <- allcomps[!allcomps %in% sigcomps]
	delta.signal <- sigex.delta(mdl,noisecomps)
	TdiffSig <- T - length(delta.signal) + 1

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	acfsignal.mat <- matrix(0,nrow=N*TdiffSig,ncol=N)
	for(i in sigcomps)
	{
		L.par[[i]] <- param[[1]][[i]]
		D.par[[i]] <- param[[2]][[i]]
		mdlType <- mdl[[2]][i]
		delta <- sigex.delta(mdl,c(noisecomps,i))
		acfsignal.mat <- acfsignal.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,
			param[[3]][[i]],N,delta,TdiffSig)		
	}
	signal.acf <- array(acfsignal.mat,dim=c(N,TdiffSig,N))

	dsig <- filter(signal,delta.signal,method="convolution",sides=1)[length(delta.signal):T,]
 	acf.sample <- acf(dsig,type="correlation",lag.max=TdiffSig,plot=FALSE)$acf
	
	alpha <- .05
	L <- TdiffSig-1
	diagnostics <- NULL
	for(i in 1:N)
	{
		auto.corr <- signal.acf[i,,i]/signal.acf[i,1,i]
		sig.corrs <- auto.corr[2:(lagall+1)]
		est.corrs <- acf.sample[2:(lagall+1),i,i]
		auto.corr <- c(rev(auto.corr),auto.corr[-1])
		all.w <- NULL

		for(h in 1:lagall) 
		{
			new.w <- sum((auto.corr[(L+1+1+h):(L+1+L)] + auto.corr[(L+1+1-h):(L+1+L-2*h)] -
				2*auto.corr[L+1+h]*auto.corr[(L+1+1):(L+1+L-h)])^2)
			all.w <- c(all.w,new.w)
		}
		ses <- qnorm(1 - alpha/2)*sqrt(abs(all.w)/TdiffSig)	
		diagnostics <- cbind(diagnostics,est.corrs + 1i*(est.corrs - sig.corrs)/ses)
	}

	return(diagnostics)
}
