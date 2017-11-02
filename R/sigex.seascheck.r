sigex.seascheck <- function(signal,period,tolerance,diffs)
{

	#################################
	#	sigex.seascheck
	#		by Tucker McElroy
	#
	#	A seasonal diagnostic, which for each series in signal
	#		gets an OLS AR fit (the signal gets differenced
 	#		(1-B)^diffs tmes), computes the roots,
	#		and checks magnitude of those close to seasonal
	#		frequencies.  Input a tolerance level for
	#		how close to seasonal frequency, and returns
	#		seasonal root magnitudes; less than 1.03 corresponds
	#		to "visual significance", and is a concern.
	#		If no root of a particular seasonal frequency
	#		is found, returns Inf.  Returns floor(period/2) values
	#		for each series 
	#
	#########################################

	x <- t(signal)
	N <- dim(x)[1]
	T <- dim(x)[2]

	mag.mat <- NULL
	for(i in 1:N) 
	{
		mags <- NULL
		diffed.x <- x[i,]
		if(diffs > 0) { diffed.x <- diff(x[i,],differences = diffs) }
		my.arfit <- ar.ols(diffed.x)
		my.roots <- polyroot(c(1,-1*my.arfit$ar[,,1]))
		polar.roots <- cbind(round(Mod(my.roots),digits=3),
			signif(2*pi/Arg(my.roots),digits=6),round(period*Arg(my.roots)/(2*pi),digits=2))
#		print(polar.roots)
		for(k in 1:(period/2))
		{
			mag <- polar.roots[abs(polar.roots[,3] - k) < tolerance,1]
			if(length(mag)==0) mag <- Inf
			if(length(mag)>1) mag <- min(mag)
			mags <- c(mags,mag)
		}
		mag.mat <- cbind(mag.mat,mags)
	}
				
	return(mag.mat)

}
