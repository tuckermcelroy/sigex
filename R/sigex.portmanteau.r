sigex.portmanteau <- function(resid,lag,nump)
{

	#################################
	#	sigex.portmanteau
	#		by Tucker McElroy
	#
	#	Computes the portmanteau statistic for resid,
	#		using lag number of acf lags, and
	#		returning chi square p value based on
	#		nump degrees of freedom (should equal number
	#		of estimated parameters in model, or length(psi))
	#
	#########################################

	x <- t(resid)
	N <- dim(x)[1]
	T <- dim(x)[2]

	dof <- lag*N^2 - nump
	if(dof <= 0) { lag <- nump + 1; dof <- 1 }
	acf.sample <- acf(t(x),type="covariance",lag.max=lag,plot=FALSE)$acf
	varinv <- solve(acf.sample[1,,])
	port <- 0
	for(h in 1:lag)
	{
		port <- port + sum(diag(acf.sample[h+1,,] %*% varinv %*% 
			t(acf.sample[h+1,,]) %*% varinv))
	}
	port <- T*port
	
	return(c(port,1-pchisq(port,df= dof)))
}



