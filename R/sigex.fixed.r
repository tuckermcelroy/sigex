sigex.fixed <- function(data,mdl,series,param,type)
{

	###############################
	#   sigex.fixed
	#	by Tucker McElroy	
	#
	#	Computes all nontrend fixed regression effects for series
	#		of a given type (string)
	#		i.e., X \hat{beta}
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	ind <- 0
	mean.mat <- NULL
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		if (k==series)
		{
			col.ind <- which(colnames(mdl[[4]][[k]])==type)
			if(length(col.ind) > 0) 
			{ 
				mean.mat <- as.matrix(reg.mat[,col.ind]) %*% 
					as.matrix(param[[4]][ind+col.ind]) 
				mean.mat <- ts(mean.mat,names=type)
			}
		}
		ind <- ind+len
	}
	
	return(mean.mat)
}
	



