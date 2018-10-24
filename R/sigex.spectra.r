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
	#	Requires: polymult, polysum, ARMAauto, VARMAauto, specFact,
	#		specFactmvar, sigex.getcycle, sigex.canonize
	#
	####################################################################

	N <- dim(as.matrix(L.par))[1]
	mdlType <- mdl[[2]][comp]
	d.delta <- length(delta)
	xi.mat <- L.par %*% diag(exp(D.par),nrow=length(D.par)) %*% t(L.par)

	#########################################################
	#   Compute VMA and VAR components for differenced models
	
	if(mdlType == "wn") 
	{ 
		cycle.order <- 0
		rho <- 0
		omega <- 0 
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		cycle.MA <- out[[2]]
		cycle.MA <- polymult(delta,cycle.MA) 
		comp.MA <- array(t(cycle.MA %x% diag(N)),c(N,N,length(cycle.MA)))
		comp.AR <- array(t(cycle.AR %x% diag(N)),c(N,N,length(cycle.AR)))
		comp.sigma <- xi.mat
	}
	if(mdlType == "canonWN")
	{
		canon.delta <- mdl[[3]][[comp]]
		canon.d <- length(canon.delta)
		freq0 <- sum(canon.delta)^2 
		freqpi <- sum(canon.delta*((-1)^(seq(0,canon.d-1))))^2
		freqmax <- max(freq0,freqpi)
		ma.poly <- polymult(rev(canon.delta),canon.delta)
		psi.acf <- c(1,rep(0,canon.d-1))
		psi.acf <- psi.acf - freqmax^{-1}*ma.poly[canon.d:(2*canon.d-1)]				
		psi.acf <- c(rev(psi.acf),psi.acf[-1])
		psi.ma <- Re(specFact(psi.acf))		
 		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
 		psi.MA <- polymult(delta,psi.ma)
		psi.AR <- 1
		comp.MA <- array(t(psi.MA %x% diag(N)),c(N,N,length(psi.MA)))
		comp.AR <- array(t(psi.AR %x% diag(N)),c(N,N,length(psi.AR)))
		comp.sigma <- psi.scale*xi.mat
	}
	if(mdlType == "AR1") 
	{
		cycle.order <- 1
		rho <- mdlPar[1]
		omega <- 0 
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		cycle.MA <- out[[2]]
		cycle.MA <- polymult(delta,cycle.MA) 
		comp.MA <- array(t(cycle.MA %x% diag(N)),c(N,N,length(cycle.MA)))
		comp.AR <- array(t(cycle.AR %x% diag(N)),c(N,N,length(cycle.AR)))
		comp.sigma <- xi.mat
	}
	if(mdlType == "MA1")
	{
		ar.poly <- 1
		ma.poly <- c(1,mdlPar[1])
		ma.poly <- polymult(delta,ma.poly)
		comp.MA <- array(t(ma.poly %x% diag(N)),c(N,N,length(ma.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- xi.mat
	}
	if(mdlType == "canonMA1")
	{
		canon.delta <- mdl[[3]][[comp]]
		psi.ma <- sigex.canonize(mdlPar[1],-1*canon.delta[-1])
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymult(delta,psi.ma)
		psi.AR <- 1
		comp.MA <- array(t(psi.MA %x% diag(N)),c(N,N,length(psi.MA)))
		comp.AR <- array(t(psi.AR %x% diag(N)),c(N,N,length(psi.AR)))
		comp.sigma <- psi.scale*xi.mat
	}
	if(mdlType %in% c("cycleBW1","cycleBW2","cycleBW3","cycleBW4","cycleBW5",
		"cycleBW6","cycleBW7","cycleBW8","cycleBW9","cycleBW10"))
	{
		if(mdlType == "cycleBW1") cycle.order <- 1
		if(mdlType == "cycleBW2") cycle.order <- 2
		if(mdlType == "cycleBW3") cycle.order <- 3
		if(mdlType == "cycleBW4") cycle.order <- 4
		if(mdlType == "cycleBW5") cycle.order <- 5
		if(mdlType == "cycleBW6") cycle.order <- 6
		if(mdlType == "cycleBW7") cycle.order <- 7
		if(mdlType == "cycleBW8") cycle.order <- 8
		if(mdlType == "cycleBW9") cycle.order <- 9
		if(mdlType == "cycleBW10") cycle.order <- 10
		rho <- mdlPar[1]
		omega <- mdlPar[2] 
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		cycle.MA <- out[[2]]
		cycle.MA <- polymult(delta,cycle.MA) 
		comp.MA <- array(t(cycle.MA %x% diag(N)),c(N,N,length(cycle.MA)))
		comp.AR <- array(t(cycle.AR %x% diag(N)),c(N,N,length(cycle.AR)))
		comp.sigma <- xi.mat
	}
	if(mdlType %in% c("canonCycleBW1","canonCycleBW2","canonCycleBW3","canonCycleBW4",
			"canonCycleBW5","canonCycleBW6","canonCycleBW7","canonCycleBW8",
			"canonCycleBW9","canonCycleBW10"))
	{
		if(mdlType == "canonCycleBW1") cycle.order <- 1
		if(mdlType == "canonCycleBW2") cycle.order <- 2
		if(mdlType == "canonCycleBW3") cycle.order <- 3
		if(mdlType == "canonCycleBW4") cycle.order <- 4
		if(mdlType == "canonCycleBW5") cycle.order <- 5
		if(mdlType == "canonCycleBW6") cycle.order <- 6
		if(mdlType == "canonCycleBW7") cycle.order <- 7
		if(mdlType == "canonCycleBW8") cycle.order <- 8
		if(mdlType == "canonCycleBW9") cycle.order <- 9
		if(mdlType == "canonCycleBW10") cycle.order <- 10
		rho <- mdlPar[1]
		omega <- mdlPar[2] 
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		cycle.MA <- out[[2]]
		psi.acf <- ARMAauto(ar = NULL,ma = cycle.MA[-1],lag.max=length(cycle.MA))
		freq0 <- ((1-2*rho*cos(pi*omega)+rho^2*cos(pi*omega)^2)/(1+rho^2-2*rho*cos(pi*omega))^2)^cycle.order
		freqpi <- ((1+2*rho*cos(pi*omega)+rho^2*cos(pi*omega)^2)/(1+rho^2+2*rho*cos(pi*omega))^2)^cycle.order
		lambda.crit1 <- (1+rho^2*cos(pi*omega)^2 - sin(pi*omega)*sqrt(sin(pi*omega)^2 + cos(pi*omega)^2*(1-rho^2)^2))/(2*rho*cos(pi*omega))
		lambda.crit2 <- (1+rho^2*cos(pi*omega)^2 + sin(pi*omega)*sqrt(sin(pi*omega)^2 + cos(pi*omega)^2*(1-rho^2)^2))/(2*rho*cos(pi*omega))
		if(abs(lambda.crit1)<=1) { lambda.crit1 <- acos(lambda.crit1) } else { lambda.crit1 <- 0 }
		if(abs(lambda.crit2)<=1) { lambda.crit2 <- acos(lambda.crit2) } else { lambda.crit2 <- 0 }
		freqcrit1 <- ((1+rho^2*cos(pi*omega)^2-2*rho*cos(pi*omega)*cos(lambda.crit1))/
			(1+4*rho^2*cos(pi*omega)^2 + rho^4 - 4*rho*(1+rho^2)*cos(pi*omega)*cos(lambda.crit1) + 2*rho^2*cos(2*lambda.crit1)))^cycle.order
		freqcrit2 <- ((1+rho^2*cos(pi*omega)^2-2*rho*cos(pi*omega)*cos(lambda.crit2))/
			(1+4*rho^2*cos(pi*omega)^2 + rho^4 - 4*rho*(1+rho^2)*cos(pi*omega)*cos(lambda.crit2) + 2*rho^2*cos(2*lambda.crit2)))^cycle.order
		freqmin <- min(freq0,freqpi,freqcrit1,freqcrit2)
		psi.acf <- c(rev(psi.acf),psi.acf[-1])
		psi.acf <- polysum(psi.acf,-freqmin*polymult(cycle.AR,rev(cycle.AR)))
		psi.ma <- Re(specFact(psi.acf))	
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymult(delta,psi.ma)
		psi.AR <- cycle.AR
		comp.MA <- array(t(psi.MA %x% diag(N)),c(N,N,length(psi.MA)))
		comp.AR <- array(t(psi.AR %x% diag(N)),c(N,N,length(psi.AR)))
		comp.sigma <- psi.scale*xi.mat
	}	
	if(mdlType %in% c("cycleBAL1","cycleBAL2","cycleBAL3","cycleBAL4","cycleBAL5",
		"cycleBAL6","cycleBAL7","cycleBAL8","cycleBAL9","cycleBAL10"))
	{
		if(mdlType == "cycleBAL1") cycle.order <- 1
		if(mdlType == "cycleBAL2") cycle.order <- 2
		if(mdlType == "cycleBAL3") cycle.order <- 3
		if(mdlType == "cycleBAL4") cycle.order <- 4
		if(mdlType == "cycleBAL5") cycle.order <- 5
		if(mdlType == "cycleBAL6") cycle.order <- 6
		if(mdlType == "cycleBAL7") cycle.order <- 7
		if(mdlType == "cycleBAL8") cycle.order <- 8
		if(mdlType == "cycleBAL9") cycle.order <- 9
		if(mdlType == "cycleBAL10") cycle.order <- 10
		rho <- mdlPar[1]
		omega <- mdlPar[2]
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		r <- seq(0,cycle.order)
		ma.acf <- sum((choose(cycle.order,r)^2)*(-rho)^(2*r))
		for(h in 1:cycle.order)
		{
			r <- seq(0,cycle.order-h)
			new.acf <- cos(h*pi*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
			ma.acf <- c(ma.acf,new.acf)
		}		
		ma.acf <- c(rev(ma.acf),ma.acf[-1])
		psi.ma <- Re(specFact(ma.acf))	
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymult(delta,psi.ma)
		psi.AR <- cycle.AR
		comp.MA <- array(t(psi.MA %x% diag(N)),c(N,N,length(psi.MA)))
		comp.AR <- array(t(psi.AR %x% diag(N)),c(N,N,length(psi.AR)))
		comp.sigma <- psi.scale*xi.mat
	}
if(mdlType %in% c("canonCycleBAL1","canonCycleBAL2","canonCycleBAL3","canonCycleBAL4",
		"canonCycleBAL5","canonCycleBAL6","canonCycleBAL7","canonCycleBAL8",
		"canonCycleBAL9","canonCycleBAL10"))
	{
		if(mdlType == "canonCycleBAL1") cycle.order <- 1
		if(mdlType == "canonCycleBAL2") cycle.order <- 2
		if(mdlType == "canonCycleBAL3") cycle.order <- 3
		if(mdlType == "canonCycleBAL4") cycle.order <- 4
		if(mdlType == "canonCycleBAL5") cycle.order <- 5
		if(mdlType == "canonCycleBAL6") cycle.order <- 6
		if(mdlType == "canonCycleBAL7") cycle.order <- 7
		if(mdlType == "canonCycleBAL8") cycle.order <- 8
		if(mdlType == "canonCycleBAL9") cycle.order <- 9
		if(mdlType == "canonCycleBAL10") cycle.order <- 10
		rho <- mdlPar[1]
		omega <- mdlPar[2]
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		r <- seq(0,cycle.order)
		ma.acf <- sum((choose(cycle.order,r)^2)*(-rho)^(2*r))
		for(h in 1:cycle.order)
		{
			r <- seq(0,cycle.order-h)
			new.acf <- cos(h*pi*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
			ma.acf <- c(ma.acf,new.acf)
		}		
		ma.acf <- c(rev(ma.acf),ma.acf[-1])
		freq0 <- (1 -2*rho*cos(pi*omega) + rho^2)^(-cycle.order)
		freqpi <- (1 +2*rho*cos(pi*omega) + rho^2)^(-cycle.order)
		freqmin <- min(freq0,freqpi)
		ma.acf <- polysum(ma.acf,-freqmin*polymult(cycle.AR,rev(cycle.AR)))
		psi.ma <- Re(specFact(ma.acf))	
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymult(delta,psi.ma)
		psi.AR <- cycle.AR
		comp.MA <- array(t(psi.MA %x% diag(N)),c(N,N,length(psi.MA)))
		comp.AR <- array(t(psi.AR %x% diag(N)),c(N,N,length(psi.AR)))
		comp.sigma <- psi.scale*xi.mat
	}
	if(mdlType == "ARMA22")
	{
		ar.poly <- c(1,-1*mdlPar[1:2])
		ma.poly <- c(1,mdlPar[3:4])
		ma.poly <- polymult(delta,ma.poly)
		comp.MA <- array(t(ma.poly %x% diag(N)),c(N,N,length(ma.poly)))
		comp.AR <- array(t(ar.poly %x% diag(N)),c(N,N,length(ar.poly)))
		comp.sigma <- xi.mat
	}
	if(mdlType == "damped")
	{
		L.zeta <- mdlPar[[1]]
		D.zeta <- mdlPar[[2]]
		phi <- mdlPar[[3]]
		phi.mat <- phi*diag(N)
		zeta.mat <- L.zeta %*% diag(exp(D.zeta),nrow=N) %*% t(L.zeta)
		varma.phi <- array(phi.mat,c(N,N,1))
		vma.acf <- array(0,c(N,2,N))
		vma.acf[,1,] <- xi.mat + zeta.mat + phi.mat %*% zeta.mat %*% t(phi.mat)
		vma.acf[,2,] <- -1*phi.mat %*% zeta.mat
		vma.theta <- specFactmvar(vma.acf)
		vma.sigma <- vma.theta[[2]]
		vma.theta <- vma.theta[[1]][,,1]
		varma.theta <- array(t(c(delta,0)) %x% diag(N) + t(c(0,delta)) %x% vma.theta,c(N,N,+1))
		comp.MA <- varma.theta
		comp.AR <- varma.phi
		comp.sigma <- vma.sigma
	}	
	
	####################################
	# compute VMA and VAR spectra

	lambda <- pi*seq(0,grid)/grid
	f.ma <- rep(1,(grid+1)) %x% comp.MA[,,1]  
	if(dim(comp.MA)[3] > 1) {
	for(i in 2:dim(comp.MA)[3])
	{
		f.ma <- f.ma + exp(-1i*lambda*(i-1)) %x% comp.MA[,,i]  
	} }
	f.ma <- array(t(f.ma),c(N,N,(grid+1)))
	f.ar <- rep(1,(grid+1)) %x% comp.AR[,,1]  
	if(dim(comp.AR)[3] > 1) {
	for(i in 2:dim(comp.AR)[3])
	{
		f.ar <- f.ar + exp(-1i*lambda*(i-1)) %x% comp.AR[,,i] 
	} }
	f.ar <- array(t(f.ar),c(N,N,(grid+1)))

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
