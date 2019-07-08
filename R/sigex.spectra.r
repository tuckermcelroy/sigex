sigex.spectra <- function(L.par,D.par,mdl,comp,mdlPar,delta,grid)
{

	##########################################################################
	#
	#	sigex.spectra
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
	#	Purpose: computes scalar part of spectrum of a differenced latent 
	#		multivariate component process
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
 	#		grid: desired number of frequencies for output spectrum
	#	Outputs:
	#		f.spec: array of dimension N x N x (grid+1), consisting of spectrum
	#			at frequencies pi*j/grid for 0 <= j <= grid
	#	Requires: polymult, polysum, polymulMat, ARMAauto, VARMAauto, 
	#		specFact, specFactmvar, sigex.getcycle, sigex.canonize
	#
	####################################################################

	N <- dim(as.matrix(L.par))[1]
	mdlType <- mdl[[2]][[comp]]
	mdlClass <- mdlType[[1]]
	mdlOrder <- mdlType[[2]]
	mdlBounds <- mdlType[[3]]
	d.delta <- length(delta)
	xi.mat <- L.par %*% diag(exp(D.par),nrow=length(D.par)) %*% t(L.par)

	#########################################################
	#   Compute VMA and VAR components for differenced models
	
	# ARMA model
	if(mdlClass == "arma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0) ar.coef <- mdlPar[1:p.order]
		if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
		ar.poly <- c(1,-1*ar.coef)
		ma.poly <- c(1,ma.coef)
		ma.poly <- polymult(delta,ma.poly)
		comp.MA <- array(t(ma.poly %x% diag(N)),c(N,N,length(ma.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- xi.mat
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
		ar.poly <- c(1,-1*ar.coef)
		canon.delta <- mdl[[3]][[comp]]
		ardiff.poly <- polymult(c(1,-1*ar.coef),canon.delta)
		ma.stab <- sigex.canonize(ma.coef,-1*ardiff.poly[-1])
		ma.scale <- ma.stab[1]^2
		ma.stab <- ma.stab/ma.stab[1]	
		madiff.stab <- polymult(delta,ma.stab)
		comp.MA <- array(t(madiff.stab %x% diag(N)),c(N,N,length(madiff.stab)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- ma.scale*xi.mat
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
		ma.poly <- polymult(ma.poly,delta)
		comp.MA <- array(t(ma.poly %x% diag(N)),c(N,N,length(ma.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- xi.mat
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
		comp.MA <- array(t(madiff.stab %x% diag(N)),c(N,N,length(madiff.stab)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- ma.scale*xi.mat
	}

	# VARMA model
	if(mdlClass == "varma")
	{
		p.order <- mdlOrder[1]
		q.order <- mdlOrder[2]
		ar.coef <- NULL
		ma.coef <- NULL
		if(p.order > 0) ar.coef <- matrix(mdlPar[,,1:p.order,drop=FALSE],nrow=N)
		if(q.order > 0) ma.coef <- matrix(mdlPar[,,(p.order+1):(p.order+q.order),drop=FALSE],nrow=N)
		ar.array <- array(cbind(diag(N),-1*ar.coef),c(N,N,p.order+1))
		ma.array <- array(cbind(diag(N),ma.coef),c(N,N,q.order+1))
		delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
		madiff.array <- polymulMat(delta.array,ma.array) 
		comp.MA <- madiff.array
		comp.AR <- ar.array
		comp.sigma <- xi.mat
	}

	# SVARMA model
	if(mdlClass == "svarma")
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
		ar.array <- array(diag(N),c(N,N,1))
		ma.array <- array(diag(N),c(N,N,1))
		ars.array <- array(diag(N),c(N,N,1))
		mas.array <- array(diag(N),c(N,N,1))
		if(p.order > 0) 
		{
			ar.coef <- mdlPar[,,1:p.order,drop=FALSE]
			ar.array <- array(cbind(diag(N),-1*matrix(ar.coef,nrow=N)),c(N,N,p.order+1))
		}
		if(q.order > 0) 
		{
			ma.coef <- mdlPar[,,(p.order+1):(p.order+q.order),drop=FALSE]
			ma.array <- array(cbind(diag(N),-1*matrix(ma.coef,nrow=N)),c(N,N,q.order+1))
		}
		if(ps.order > 0) 
		{
			ars.coef <- matrix(mdlPar[,,(p.order+q.order+1):(p.order+q.order+ps.order),drop=FALSE],c(N,N,ps.order))
			ars.coef <- array(t(stretch) %x% ars.coef,c(N,N,s.period*ps.order)) 
			ars.array <- array(cbind(diag(N),-1*matrix(ars.coef,nrow=N)),c(N,N,s.period*ps.order+1))
		}
		if(qs.order > 0)
		{
			mas.coef <- matrix(mdlPar[,,(p.order+q.order+ps.order+1):(p.order+q.order+ps.order+qs.order),drop=FALSE],c(N,N,qs.order))
			mas.coef <- array(t(stretch) %x% mas.coef,c(N,N,s.period*qs.order))
			mas.array <- array(cbind(diag(N),-1*matrix(mas.coef,nrow=N)),c(N,N,s.period*qs.order+1))
		}
		ar.poly <- polymulMat(ar.array,ars.array) 
		ma.poly <- polymulMat(ma.array,mas.array) 
		delta.array <- array(t(delta) %x% diag(N),c(N,N,d.delta))
		madiff.array <- polymulMat(delta.array,ma.poly) 
		comp.MA <- madiff.array
		comp.AR <- ar.poly
		comp.sigma <- xi.mat
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
		madiff.poly <- polymult(delta,ma.poly) 
		comp.MA <- array(t(madiff.poly %x% diag(N)),c(N,N,length(madiff.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- xi.mat
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
		comp.MA <- array(t(madiff.stab %x% diag(N)),c(N,N,length(madiff.stab)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- ma.scale*xi.mat
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
		ma.poly <- polymult(delta,ma.poly)
		comp.MA <- array(t(ma.poly %x% diag(N)),c(N,N,length(ma.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- ma.scale*xi.mat
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
		comp.MA <- array(t(madiff.stab %x% diag(N)),c(N,N,length(madiff.stab)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- ma.scale*xi.mat
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
		ma.poly <- delta.poly
		comp.MA <- array(t(ma.poly %x% diag(N)),c(N,N,length(ma.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- xi.mat
	}	

	#############################
	# compute VMA and VAR spectra

	lambda <- pi*seq(0,grid)/grid
	f.ma <- t(rep(1,(grid+1))) %x% comp.MA[,,1]  
	if(dim(comp.MA)[3] > 1) {
	for(i in 2:dim(comp.MA)[3])
	{
		f.ma <- f.ma + t(exp(-1i*lambda*(i-1))) %x% comp.MA[,,i]  
	} }
	f.ma <- array(f.ma,c(N,N,(grid+1)))
	f.ar <- t(rep(1,(grid+1))) %x% comp.AR[,,1]  
	if(dim(comp.AR)[3] > 1) {
	for(i in 2:dim(comp.AR)[3])
	{
		f.ar <- f.ar + t(exp(-1i*lambda*(i-1))) %x% comp.AR[,,i] 
	} }
	f.ar <- array(f.ar,c(N,N,(grid+1)))

	### compute spectrum
	f.wold <- array(t(rep(1,(grid+1)) %x% diag(N)),c(N,N,(grid+1)))
	f.spec <- f.wold
	for(k in 1:(grid+1))
	{
		f.wold[,,k] <- solve(f.ar[,,k]) %*% f.ma[,,k]
	}
	for(k in 1:(grid+1))
	{
		f.spec[,,k] <- f.wold[,,k] %*% comp.sigma %*% Conj(t(f.wold[,,k]))
	}

	return(f.spec)
}
