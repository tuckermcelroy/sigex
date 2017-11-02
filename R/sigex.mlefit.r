sigex.mlefit <- function(data,param,flag,mdl,method,hess=TRUE)
{

	###############################
	#   SigEx mlefit
	#	by Tucker McElroy	
	#
	#	Fits mdl to data x, with initialization of param
	#		param must be in format yielded by sigex.default
	#		To fix values, add 1i to it, as a flag
	#	   Method is "bfgs" for BFGS, "sann" for simulated annealing
	#		hess is a Boolean flag; if true, for BFGS it runs
	#		another round of BFGS to get Hessian matrix
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

	# set thresholding to prevent crash (irrelevant for SANN)
	lower.bound <- rep(-30,length(psi.est))
	upper.bound <- rep(10,length(psi.est))

	# initial attempt to fit
	if(method=="bfgs") {
	mle <- try(nlminb(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data=data,flag=flag,
		lower=lower.bound,upper=upper.bound,
		control=list(iter.max=50,eval.max=200)),TRUE) 
	}
	if(method=="sann") {
	mle <- optim(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data=data,
		flag=flag,method="SANN",control=list(maxit=5000)) }

	# if the fit was successful, we can report results;
	#	also, we can try to get the hessian - in this case,
	#	we overwrite previous results with L-BFGS-B method
	#	(if not desired, set hess = FALSE)
	if(!inherits(mle, "try-error")) {
		if(hess) {
			psi.est <- mle$par
			mle.hess <- try(optim(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,
				data=data,flag=flag,
				lower=lower.bound,upper=upper.bound,
				method="L-BFGS-B",hessian=TRUE,
				control=list(maxit=100))) 
			if(!inherits(mle.hess, "try-error")) { mle <- mle.hess }
		}
		psi.est <- mle$par
		psi <- flag
		psi[flag==1] <- psi.est
		psi[flag==0] <- psi.fix
		par.est <- sigex.psi2par(psi,mdl,data) } else mle <- Inf

	return(list(mle,par.est))
}
