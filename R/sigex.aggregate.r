sigex.aggregate <- function(data,filter,mdl,param,beta)
{

	###############################
	#   sigex.aggregate
	#	by Tucker McElroy	
	#
	#	Computes signal extraction estimates with two standard errors
	#		    for some aggregate of component signals
	# 		filter is list object output of sigex.signal
	#		param must be in format yielded by sigex.default
	#		beta is N vector of linear combinations of signal components
	#     Output is list of object with three length T vectors: the estimate,
	#		followed by upper and lower confidence intervals
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# subtract regression effects
	ind <- 0
	data.diff <- data
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		data.diff[,k] <- data.diff[,k] - reg.mat %*% 
			as.matrix(param[[4]][(ind+1):(ind+len)])
		ind <- ind+len
	}
	xvec <- matrix(t(data.diff),nrow=N*T,ncol=1)

#	trendComp <- sigex.whichtrend(mdl)
#	allcomps <- seq(1,length(mdl[[3]]))
#	dOrd <- length(mdl[[3]][[trendComp]])-1
#	deltaNontrend <- sigex.delta(mdl,trendComp)
#	period <- sum(deltaNontrend)

#	xvec <- matrix(t(x),nrow=N*T,ncol=1)
#	if(dOrd > 0)
#	{
#		meanV <- param[[5]]
#		nu <- meanV/(factorial(dOrd)*period)
#		tau <- seq(1,T)^dOrd
#		meanFcn <- matrix(nu,ncol=1) %x% matrix(tau,ncol=1)	
#		xvec <- xvec - meanFcn
#	}

   	extract <- filter[[1]] %*% xvec
	extract <- t(matrix(extract,nrow=N,ncol=T))
	extract <- extract %*% beta

	mse <- (t(beta) %x% diag(T)) %*% filter[[2]] %*% t(t(beta) %x% diag(T))
	mse <- diag(mse)
	upp <- extract + 2*sqrt(mse)
	low <- extract - 2*sqrt(mse)
	 
	return(list(extract,upp,low))
}
