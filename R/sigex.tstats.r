sigex.tstats <- function(mdl,psi,hess)
{

	############################################
	#
	#	SigEx Tstats
	#		Tucker McElroy
	#
	#	Computes standard errors for parameter estimates,
	#		and returns tstats
	#	NOTES:
	#     standard error has no division by T, because mvar.lik
	#	is T times the scale of the Whittle lik; mult by 2 because 
	#	mvar.lik is -2* log lik
	#
	#####################################

	se <- rep(0,sum(Im(psi)))
	if(min(eigen(hess)$value) > 0) se <- sqrt(2*diag(solve(hess)))
	
	tstats <- Re(psi)
	tstats[Im(psi)==1] <- tstats[Im(psi)==1]/se
	tstats[Im(psi)==0] <- sign(tstats[Im(psi)==0])*Inf
	return(tstats)
}


