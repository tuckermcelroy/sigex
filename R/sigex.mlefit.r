sigex.mlefit <- function(data.ts,param,constraint,mdl,method,thresh=Inf,hess=TRUE,whittle=FALSE,debug=FALSE)
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
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter manifold)  
	#	Notes: handles missing values in data.ts, which are indicated by 1i.
	#	Inputs:
	#		data.ts: a T x N matrix ts object; any missing values 
	#			must be encoded with 1i in that entry
	#		param: see background; this is an initial specification to 
	#			start the nonlinear optimization routines
  #		constraint: matrix of the form [Q , C], with C (constraint.mat)
  #     the matrix of constraints and Q (constraint.vec) the vector
  #     of constraint constants, such that C psi = Q. 
  #     Use NULL if there are no constraints
	#		mdl: the specified sigex model, a list object
	#	  method: "bfgs" for BFGS, "sann" for simulated annealing, 
  #     "cg" for conjugate gradient
  #   thresh: pre-parameters theta satisfy |theta|< thresh; 
  #     set thresh = Inf if no thresholding is desired
	#		hess: a Boolean flag; if true, for BFGS it runs
	#			another round of BFGS to get Hessian matrix
	#		whittle: a Boolean flag; if true, uses Whittle likelihood instead of
	#			default Gaussian likelihood	
	#		debug: a Boolean flag; if true, sets the DEBUGGING mode, which prints
	#			psi.last and psi.now to the global field.  Also lik will be printed.
	#			psi.now gives the parameter that crashed the likelihood (if it crashed),
	#			psi.last gives the last good parameter that did not crash the lik.
	#	Outputs:
	#		list with mle and par.est
	#		mle: an object of type outputted by optim
	#		par.est: type param, with the estimated parameters filled in
	#	Requires: sigex.par2psi, sigex.psi2par, sigex.lik, sigex.whittle, 
  #   sigex.eta2psi, sigex.psi2eta
	#
	####################################################################
  
  fix.lik <- function(eta,constraint,mdl,data.ts,whittle,debug)
  {
    psi <- sigex.eta2psi(eta,constraint)
    if(debug) psi.now <<- psi
    if(whittle)
    {
      out <- sigex.whittle(psi,mdl,data.ts)
    } else
    {
      out <- sigex.lik(psi,mdl,data.ts,debug)
    }
    if(debug) psi.last <<- psi
    return(out)
  }
  
  x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	par.est <- NULL
	psi <- sigex.par2psi(param,mdl)
	# check: should be zero if constraint holds for initial value
#	if(length(constraint)>0)
#	{
#	  check <- sum((constraint[,-1] %*% psi - constraint[,1])^2)
#	  print(check)
#	}
  nueta <- sigex.psi2eta(psi,constraint)
  nu <- nueta[[1]]
  eta <- nueta[[2]]

	# set thresholding to prevent crash (irrelevant for SANN)
	lower.bound <- rep(-thresh,length(eta))
	upper.bound <- rep(thresh,length(eta))
	
	# initial attempt to fit
	if(method=="bfgs") {
	mle <- try(nlminb(eta,fix.lik,constraint=constraint,mdl=mdl,data.ts=data.ts,
	          whittle=whittle,debug=debug,lower=lower.bound,upper=upper.bound,
		        control=list(iter.max=50,eval.max=200)),TRUE) 
	}
	if(method=="sann") {
	mle <- optim(eta,fix.lik,constraint=constraint,mdl=mdl,data.ts=data.ts,
	             whittle=whittle,debug=debug,method="SANN",control=list(maxit=5000)) }
	if(method=="cg") {
	mle <- optim(eta,fix.lik,constraint=constraint,mdl=mdl,data.ts=data.ts,
	             whittle=whittle,debug=debug,method="CG",control=list(maxit=50,trace=10)) }

	# if the fit was successful, we can report results;
	#	also, we can try to get the hessian - in this case,
	#	we overwrite previous results with L-BFGS-B method
	#	(if not desired, set hess = FALSE)
	if(!inherits(mle, "try-error")) {
		if(hess) {
			eta.est <- mle$par
			mle.hess <- try(optim(eta,fix.lik,constraint=constraint,mdl=mdl,data.ts=data.ts,
			                      whittle=whittle,debug=debug,lower=lower.bound,upper=upper.bound,
				                    method="L-BFGS-B",hessian=TRUE,control=list(maxit=100))) 
			if(!inherits(mle.hess, "try-error")) { mle <- mle.hess }
		}
		eta.est <- mle$par
		psi.est <- sigex.eta2psi(eta.est,constraint)
		par.est <- sigex.psi2par(psi.est,mdl,data.ts) } else mle <- Inf

	return(list(mle,par.est))
}
