sigex.mlehess <- function(data,param,flag,mdl)
{

	###############################
	#   SigEx mlehess
	#	by Tucker McElroy	
	#
	#	Fits mdl to data, with initialization of param
	#		param must be in format yielded by sigex.default
	#		To fix values, add 1i to it, as a flag
	#	   Method is  BFGS, using optim and yields numerical Hessian
	#		utilize based on initialization from a previous fit
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	par.est <- NULL
	psi <- sigex.par2psi(param,flag,mdl)
	psi <- Re(psi)
	psi.est <- psi[flag==1]
	psi.fix <- psi[flag==0]	
	fix.lik <- function(psi.est,psi.fix,mdl,data,flag)
	{
		psi.full <- flag
		psi.full[flag==1] <- psi.est
		psi.full[flag==0] <- psi.fix
		out <- sigex.lik(psi.full,mdl,data)
		return(out)
	}
#	hess <- optim(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data=data,
#		flag=flag,method="BFGS",hessian=TRUE,control=list(maxit=1)) 
	hess <- optimHess(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data=data,
		flag=flag)

	return(hess)
}
