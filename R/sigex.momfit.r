sigex.momfit <- function(data.ts,param,mdl)
{

	##########################################################################
	#
	#	sigex.momfit
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
	#	Purpose: computes initial parameter estimates by method of moments
	#	Background:
	#		param is the name for the model parameters entered into
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Notes: does not handle missing values in data.ts.
	#		ARMA model parameters are not estimated, but taken as given
	#		Does not work with VARMA/SVARMA specification
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background; this is an initial specification for
	#			the covariance matrix parameters, taking all ARMA
	#			parameters as fixed.
	#		mdl: the specified sigex model, a list object
	#	Outputs:
	#		par.new: type param, with the estimated covariance parameters filled in
	#	Requires: sigex.delta, sigex.getcycle, polymult, polysum, specFact,
	#		ARMAauto, getGCD, sigex.canonize
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]
	par.new <- param

	fulldiff <- sigex.delta(mdl,0)
	data.diff <- filter(t(x),fulldiff,method="convolution",
			sides=1)[length(fulldiff):T,,drop=FALSE]
	Tdiff <- dim(data.diff)[1]

	# get OLS estimates of regressors
	betas.ols <- NULL
	mu.ols <- NULL
	for(k in 1:N)
	{
		reg <- mdl[[4]][[k]]
		reg.diff <- as.matrix(filter(reg,fulldiff,method="convolution",
			sides=1)[length(fulldiff):T,])
		reg.mat <- t(reg.diff) %*% reg.diff
		reg.y <- t(reg.diff) %*% data.diff[,k]
		beta.ols <- solve(reg.mat) %*% reg.y
		# significance thresholding
		resid <- data.diff[,k] - reg.diff %*% beta.ols
		error.cov <- solve(reg.mat) * sum(resid^2)/length(resid)
		beta.ols[beta.ols < 2*sqrt(diag(error.cov))] <- 0
		betas.ols <- rbind(betas.ols,beta.ols)
		mu.ols <- cbind(mu.ols,reg.diff%*%beta.ols)
	}

	mu <- mu.ols
	data.diff <- data.diff - mu
	data.acf <- acf(data.diff,plot=FALSE,lag.max=T-1,type="covariance")$acf

	ma.pols <- NULL
	ar.pols <- NULL
	for(i in 1:length(mdl[[3]]))
	{
		delta.poly <- sigex.delta(mdl,i)
		mdlType <- mdl[[2]][[i]]
		mdlClass <- mdlType[[1]]
		mdlOrder <- mdlType[[2]]
		mdlBounds <- mdlType[[3]]
		mdlPar <- param[[3]][[i]]

		# ARMA model
		if(mdlClass == "arma")
		{
			p.order <- mdlOrder[1]
			q.order <- mdlOrder[2]
			ar.coef <- NULL
			ma.coef <- NULL
			if(p.order > 0) ar.coef <- mdlPar[1:p.order]
			if(q.order > 0) ma.coef <- mdlPar[(p.order+1):(p.order+q.order)]
			ma.poly <- polymult(c(1,ma.coef),delta.poly)
			ar.poly <- c(1,-1*ar.coef)
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
			canon.delta <- mdl[[3]][[i]]
			ardiff.poly <- polymult(c(1,-1*ar.coef),canon.delta)
			ma.stab <- sigex.canonize(ma.coef,-1*ardiff.poly[-1])
			ma.poly <- polymult(delta.poly,ma.stab)
			ar.poly <- c(1,-1*ar.coef)
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
			ma.poly <- polymult(ma.poly,delta.poly)
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
			ma.poly <- polymult(delta.poly,ma.stab)
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
			ma.poly <- polymult(delta.poly,ma.poly)
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
			canon.delta <- mdl[[3]][[i]]
			ardiff.poly <- polymult(ar.poly,canon.delta)
			ma.stab <- sigex.canonize(ma.poly[-1],-1*ardiff.poly[-1])
			ma.poly <- polymult(delta.poly,ma.stab)
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
			ma.poly <- polymult(delta.poly,ma.poly)
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
			ma.scale <- abs(ma.poly[1])
			ma.coef <- ma.poly[-1]/ma.poly[1]
			canon.delta <- mdl[[3]][[i]]
			ardiff.poly <- polymult(ar.poly,canon.delta)
			ma.stab <- sigex.canonize(ma.coef,-1*ardiff.poly[-1])
			ma.stab <- ma.scale*ma.stab
			ma.poly <- polymult(delta,ma.stab)
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
		}

		ma.null <- NULL
		if(length(ma.poly)==1) ma.null <- 0
		ma.pols[[length(ma.pols)+1]] <- c(ma.null,ma.poly)
		if(length(ma.poly)==1) ma.pols[[length(ma.pols)]] <- ma.pols[[length(ma.pols)]][-1]
		ar.null <- NULL
		if(length(ar.poly)==1) ar.null <- 0
		ar.pols[[length(ar.pols)+1]] <- c(ar.null,ar.poly)
		if(length(ar.poly)==1) ar.pols[[length(ar.pols)]] <- ar.pols[[length(ar.pols)]][-1]
	}

	est.acf <- array(0,c(N,length(mdl[[3]]),N))
	G.mat <- matrix(0,length(mdl[[3]]),length(mdl[[3]]))
	for(i in 1:length(mdl[[3]]))
	{
		for(j in 1:length(mdl[[3]]))
		{
			ma.prod <- polymult(ma.pols[[i]],ma.pols[[j]])
			ma.scale <- ma.prod[1]^2
			ma.prod <- ma.prod/ma.prod[1]
			ar.prod <- polymult(ar.pols[[i]],ar.pols[[j]])
			G.mat[i,j] <- ARMAauto(ar = -1*ar.prod[-1],
				ma = ma.prod[-1],lag.max=1)[1]*ma.scale
		}
		g.acf <- ARMAauto(ar = -1*ar.pols[[i]][-1],ma = ma.pols[[i]][-1],lag.max=Tdiff-1)
		new.acf <- g.acf[1]*data.acf[1,,]
		for(k in 2:Tdiff)
		{
			new.acf <- new.acf + g.acf[k]*(data.acf[k,,] + t(data.acf[k,,]))
		}
		est.acf[,i,] <- new.acf
	}
	G.mat.inv <- solve(G.mat)
	par.est.mat <- G.mat.inv %x% diag(N) %*% matrix(est.acf,c(length(mdl[[3]])*N,N))

	for(i in 1:length(mdl[[3]]))
	{
		temp.mat <- par.est.mat[((i-1)*N+1):(i*N),]
		if(N==1) new.mat <- as.matrix(pmax(0,temp.mat))	else {
			eig.mat <- eigen(temp.mat)
			new.mat <- eig.mat$vectors %*% diag(pmax(0,eig.mat$values)) %*% t(eig.mat$vectors)
		}
		gcd.est <- getGCD(new.mat,N)
		par.new[[1]][[i]] <- gcd.est[[1]]
		par.new[[2]][[i]] <- log(gcd.est[[2]])
	}
	par.new[[4]] <- betas.ols

	return(par.new)
}

