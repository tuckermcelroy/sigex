sigex.acf <- function(L.par,D.par,mdl,comp,mdlPar,delta,maxlag)
{
	
	##########################################################################
	#
	#	sigex.acf
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
	#	Purpose: compute the autocovariance function of a differenced latent component
	#	Background:	
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		We have a model for each w_t process, and can compute its
	#		autocovariance function (acf), and denote its autocovariance
	#		generating function (acgf) via gamma_w (B).
	#			Sometimes we may over-difference,
	#		which means applying a differencing polynomial eta(B) that
	#		contains delta(B) as a factor: eta(B) = delta(B)*nu(B).
	#		Then  eta(B) y_t = nu(B) w_t, and the corresponding 
	#		acgf is   nu(B) * nu(B^{-1}) * gamma_w (B). 
	#	Notes: this function computes the over-differenced acgf, 
	#		it is presumed that the given eta(B) contains the needed delta(B)
	#		for that particular component.
	#	Inputs:
	#		L.par: unit lower triangular matrix in GCD of the component's	
	#			white noise covariance matrix.  (Cf. sigex.param2gcd background)
	#		D.par: vector of logged entries of diagonal matrix in GCD
	#			of the component's white noise covariance matrix.
	#			(Cf. sigex.param2gcd background)
	#		mdl: the specified sigex model, a list object
	#		comp: index of the latent component
	#		mdlPar: see background to sigex.par2zeta.  This is the portion of param
	#			corresponding to mdl[[2]], cited as param[[3]]
	#		delta: differencing polynomial (corresponds to eta(B) in Background)
	#			written in format c(delta0,delta1,...,deltad)
 	#		maxlag: number of autocovariances required
	#	Outputs:
	#		x.acf: matrix of dimension N x N*maxlag, consisting of autocovariance
	#			matrices stacked horizontally, i.e.
	#			x.acf = [ gamma(0), gamma(1), ..., gamma(maxlag-1)]
	#	Requires: polymult, polysum, ARMAauto, VARMAauto, specFact,
	#		specFactmvar, sigex.getcycle, sigex.canonize
	#
	####################################################################

	mdlType <- mdl[[2]][[comp]]
	mdlClass <- mdlType[[1]]
	mdlOrder <- mdlType[[2]]
	mdlBounds <- mdlType[[3]]
	d.delta <- length(delta)
	xi.mat <- L.par %*% diag(exp(D.par),nrow=length(D.par)) %*% t(L.par)
	
	##################################
	## get acf of stationary component

	# ARMA model
	if(mdlClass == "arma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		psi.acf <- ARMAauto(ar = ar.coef, ma = polymult(c(1,ma.coef),delta)[-1], 
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}

	# Stabilized ARMA model
	if(mdlClass == "arma.stab")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(c(1,-1*ar.coef),canon.delta)
		ma.stab <- sigex.canonize(ma.coef,-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]	
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = ar.coef, ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]	
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# SARMA model
	if(mdlClass == "sarma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		stretch <- c(rep(0,s.period-1),1)
		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		if(ps.order > 0) 
		{
			ars.coef <- mdlPar[(p.order+q.order+1):(p.order+q.order+ps.order)]
			ars.coef <- ars.coef %x% stretch
		}
		if(qs.order > 0)
		{
			mas.coef <- mdlPar[(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order)]
			mas.coef <- mas.coef %x% stretch
		}
		ar.poly <- polymult(c(1,-1*ar.coef),c(1,-1*ars.coef))
		ma.poly <- polymult(c(1,-1*ma.coef),c(1,-1*mas.coef))
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = polymult(ma.poly,delta)[-1], 
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}

	# Stabilized SARMA model
	if(mdlClass == "sarma.stab")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ps.order <- mdlOrder[3]
		qs.order <- mdlOrder[4]
		s.period <- mdlOrder[5]
		stretch <- c(rep(0,s.period-1),1)
		ar.coef <- NULL
		ma.coef <- NULL
		ars.coef <- NULL
		mas.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		if(ps.order > 0) 
		{
			ars.coef <- mdlPar[(p.order+q.order+1):(p.order+q.order+ps.order)]
			ars.coef <- ars.coef %x% stretch
		}
		if(qs.order > 0)
		{
			mas.coef <- mdlPar[(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order)]
			mas.coef <- mas.coef %x% stretch
		}
		ar.poly <- polymult(c(1,-1*ar.coef),c(1,-1*ars.coef))
		ma.poly <- polymult(c(1,-1*ma.coef),c(1,-1*mas.coef))
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(ar.poly,canon.delta)
		ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]	
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]	
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# Butterworth cycle
	if(mdlClass == "bw")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2] 
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		ma.poly <- out[[2]]
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1],ma = polymult(ma.poly,delta)[-1],
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}

	# Stabilized Butterworth cycle
	if(mdlClass == "bw.stab")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2] 
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		ma.poly <- out[[2]]
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(ar.poly,canon.delta)
		ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]	
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]	
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# Balanced cycle
	if(mdlClass == "bal")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2] 
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		r <- seq(0,cycle.order)
		ma.acf <- sum((choose(cycle.order,r)^2)*(-rho)^(2*r))
		for(h in 1:cycle.order)
		{
			r <- seq(0,cycle.order-h)
			new.acf <- cos(h*pi*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
			ma.acf <- c(ma.acf,new.acf)
		}		
		ma.acf <- c(rev(ma.acf),ma.acf[-1]) + 1e-10
		ma.poly <- Re(specFact(ma.acf))	
		ma.scale <- ma.poly[1]^2
		ma.poly <- ma.poly/ma.poly[1]	
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1],ma = polymult(ma.poly,delta)[-1],
			lag.max=maxlag)[1:maxlag]
	 	psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}

	# Stabilized Balanced cycle
	if(mdlClass == "bal.stab")
	{
		cycle.order <- mdlOrder[1]
		rho <- mdlPar[1]
		omega <- mdlPar[2] 
		out <- sigex.getcycle(cycle.order,rho,omega)
		ar.poly <- out[[1]]
		r <- seq(0,cycle.order)
		ma.acf <- sum((choose(cycle.order,r)^2)*(-rho)^(2*r))
		for(h in 1:cycle.order)
		{
			r <- seq(0,cycle.order-h)
			new.acf <- cos(h*pi*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
			ma.acf <- c(ma.acf,new.acf)
		}		
		ma.acf <- c(rev(ma.acf),ma.acf[-1]) + 1e-10
		ma.poly <- Re(specFact(ma.acf))	
		ma.scale <- ma.poly[1]^2
		ma.poly <- ma.poly/ma.poly[1]	
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(ar.poly,canon.delta)
		ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
		ma.scale <- ma.scale*ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]	
		madiff.stab <- polymult(delta,ma.stab)
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = madiff.stab[-1],lag.max=maxlag)[1:maxlag]	
		psi.acf <- psi.acf*ma.scale
		x.acf <- psi.acf %x% xi.mat
	}
	
	# Damped Trend model
	if(mdlClass == "damped")
	{
		p.order <- mdlOrder[1]
		ar.coef <- mdlPar[1]
		ar.poly <- 1
		for(k in 1:p.order)
		{
			ar.poly <- polymult(ar.poly,c(1,-1*ar.coef))
		}
		psi.acf <- ARMAauto(ar = -1*ar.poly[-1], ma = delta[-1],lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}	

	return(x.acf)
}
