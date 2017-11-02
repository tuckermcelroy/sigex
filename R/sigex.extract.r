sigex.extract <- function(data,filter,mdl,param)
{

	###############################
	#   sigex.extract
	#	by Tucker McElroy	
	#
	#	Computes signal extraction estimates with two standard errors
	# 		filter is list object output of sigex.signal
	#		param must be in format yielded by sigex.default
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
	 
   	extract <- filter[[1]] %*% xvec
	extract <- t(matrix(extract,nrow=N,ncol=T))

	mse <- t(matrix(diag(filter[[2]]),N,T))
	upp <- extract + 2*sqrt(mse)
	low <- extract - 2*sqrt(mse)
	 
	return(list(extract,upp,low))
}
