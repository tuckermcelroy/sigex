sigex.mlefit <- function(data.ts,param,flag,mdl,method,hess=TRUE,whittle=FALSE)
{

	##########################################################################
	#
	#	sigex.mlefit
	# 	    Copyright (C) 2017  Tucker McElroy
	#
	#    This program is free software: you can redistribute it and/or modify
	#    it under the terms of the GNU General Public License as published by
	#    the Free Software Foundation, either version 3 of the License, or
	#    (at your option) any later version.
	#
	#    This program is distributed in the hope that it will be useful,
	#    but WITHOUT ANY WARRANTY; without even the implied warranty of
	#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#    GNU General Public License for more details.
	#
	#    You should have received a copy of the GNU General Public License
	#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	#
	############################################################################

	################# Documentation #####################################
	#
	#	Purpose: fit model to the data
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter
	#		manifold) together with imaginary component flagging 
	#		whether the hyper-parameter is fixed for purposes of estimation.
	#	Notes: handles missing values in data.ts, which are indicated by 1i.
	#		To fix values of psi, set the imaginary part to zero 
	#		by modifying flag (change 1 to 0)
	#	Inputs:
	#		data.ts: a T x N matrix ts object; any missing values 
	#			must be encoded with 1i in that entry
	#		param: see background; this is an initial specification to 
	#			start the nonlinear optimization routines
	#		flag: string of numbers, of length same as psi,
	#			with a 1 denoting that the corresponding hyper-parameter is to
	#			be estimated, and 0 if it is fixe.
	#		mdl: the specified sigex model, a list object
	#	   	method: "bfgs" for BFGS, "sann" for simulated annealing, "cg" for conjugate gradient
	#		hess: a Boolean flag; if true, for BFGS it runs
	#			another round of BFGS to get Hessian matrix
	#		whittle: a Boolean flag; if true, uses Whittle likelihood instead of
	#			default Gaussian likelihood	
	#	Outputs:
	#		list with mle and par.est
	#		mle: an object of type outputted by optim
	#		par.est: type param, with the estimated parameters filled in
	#	Requires: sigex.par2psi, sigex.psi2par, sigex.lik, sigex.whittle
	#
	####################################################################
  
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	par.est <- NULL
	psi <- sigex.par2psi(param,flag,mdl)
	psi <- Re(psi)
	psi.est <- psi[flag==1]
	psi.fix <- psi[flag==0]	
	fix.lik <- function(psi.est,psi.fix,mdl,data.ts,flag)
	{
		psi.full <- flag
		psi.full[flag==1] <- psi.est
		psi.full[flag==0] <- psi.fix
		if(whittle)
		{
			out <- sigex.whittle(psi.full,mdl,data.ts)
		} else
		{
			out <- sigex.lik(psi.full,mdl,data.ts)
		}
		psi.last <<- psi.full
		return(out)
	}

	# set thresholding to prevent crash (irrelevant for SANN)
	lower.bound <- rep(-30,length(psi.est))
	upper.bound <- rep(10,length(psi.est))

	# initial attempt to fit
	if(method=="bfgs") {
	mle <- try(nlminb(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data.ts=data.ts,flag=flag,
		lower=lower.bound,upper=upper.bound,
		control=list(iter.max=50,eval.max=200)),TRUE) 
	}
	if(method=="sann") {
	mle <- optim(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data.ts=data.ts,
		flag=flag,method="SANN",control=list(maxit=5000)) }
	if(method=="cg") {
	mle <- optim(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,data.ts=data.ts,
		flag=flag,method="CG",control=list(maxit=50,trace=10)) }

	# if the fit was successful, we can report results;
	#	also, we can try to get the hessian - in this case,
	#	we overwrite previous results with L-BFGS-B method
	#	(if not desired, set hess = FALSE)
	if(!inherits(mle, "try-error")) {
		if(hess) {
			psi.est <- mle$par
			mle.hess <- try(optim(psi.est,fix.lik,psi.fix=psi.fix,mdl=mdl,
				data.ts=data.ts,flag=flag,
				lower=lower.bound,upper=upper.bound,
				method="L-BFGS-B",hessian=TRUE,
				control=list(maxit=100))) 
			if(!inherits(mle.hess, "try-error")) { mle <- mle.hess }
		}
		psi.est <- mle$par
		psi <- flag
		psi[flag==1] <- psi.est
		psi[flag==0] <- psi.fix
		par.est <- sigex.psi2par(psi,mdl,data.ts) } else mle <- Inf

	return(list(mle,par.est))
}
