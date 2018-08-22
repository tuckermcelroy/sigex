sigex.momfit <- function(data.ts,param,mdl)
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
	#	Purpose: computes initial parameter estimates by method of moments
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter
	#		manifold) together with imaginary component flagging 
	#		whether the hyper-parameter is fixed for purposes of estimation.
	#	Notes: does not handle missing values in data.ts.
	#		ARMA model parameters are not estimated, but taken as given
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background; this is an initial specification for
	#			the covariance matrix parameters, taking all ARMA
	#			parameters as fixed.
	#		mdl: the specified sigex model, a list object
	#	Outputs:
	#		par.new: type param, with the estimated covariance parameters filled in
	#	Requires: sigex.delta, sigex.getcycle, polymult, polysum, specFact,
	#		ARMA2acf, getGCD, sigex.canonize		
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
		ma.pol <- sigex.delta(mdl,i)
		ar.pol <- 1
		mdlType <- mdl[[2]][i]
	 
		if(mdlType == "canonWN")
		{
			canon.delta <- mdl[[3]][[i]]
			canon.d <- length(canon.delta)
			freq0 <- sum(canon.delta)^2 
			freqpi <- sum(canon.delta*((-1)^(seq(0,canon.d-1))))^2
			freqmax <- max(freq0,freqpi)
			ma.poly <- polymult(rev(canon.delta),canon.delta)
			psi.acf <- c(1,rep(0,canon.d-1))
			psi.acf <- psi.acf - freqmax^{-1}*ma.poly[canon.d:(2*canon.d-1)]				
			psi.acf <- c(rev(psi.acf),psi.acf[-1])
			psi.ma <- Re(specFact(psi.acf))		
 			ma.pol <- polymult(psi.ma,ma.pol)
		}	
	 	if(mdlType == "AR1") 
		{ 
			phi <- param[[3]][[i]][1]
			ar.pol <- polymult(c(1,-1*phi),ar.pol)
		}
		if(mdlType == "MA1")
		{
			theta <- param[[3]][[i]][1]
			ma.pol <- polymult(c(1,theta),ma.pol)
		}
		if(mdlType == "canonMA1")
		{
			canon.delta <- mdl[[3]][[i]]
			theta <- param[[3]][[i]][1]
			psi.ma <- sigex.canonize(theta,-1*canon.delta[-1])
			ma.pol <- polymult(psi.ma,ma.pol)
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
			rho <- param[[3]][[i]][1]
			omega <- param[[3]][[i]][2]
			out <- sigex.getcycle(cycle.order,rho,omega)
			cycle.AR <- out[[1]]
			cycle.MA <- out[[2]]
			ma.pol <- polymult(cycle.MA,ma.pol)
			ar.pol <- polymult(cycle.AR,ar.pol)
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
			rho <- param[[3]][[i]][1]
			omega <- param[[3]][[i]][2]
			out <- sigex.getcycle(cycle.order,rho,omega)
			cycle.AR <- out[[1]]
			cycle.MA <- out[[2]]
			psi.acf <- ARMA2acf(ar = NULL,ma = cycle.MA[-1],lag.max=length(cycle.MA))
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
 			ma.pol <- polymult(psi.ma,ma.pol)
			ar.pol <- polymult(cycle.AR,ar.pol)
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
				new.acf <- cos(h*omega) * sum(choose(cycle.order,r+h)*choose(cycle.order,r)*(-rho)^(2*r+h))
				ma.acf <- c(ma.acf,new.acf)
			}		
			ma.acf <- c(rev(ma.acf),ma.acf[-1])
			psi.ma <- Re(specFact(ma.acf))	
			ma.pol <- polymult(psi.ma,ma.pol)
			ar.pol <- polymult(cycle.AR,ar.pol)
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
			ma.pol <- polymult(psi.ma,ma.pol)
			ar.pol <- polymult(cycle.AR,ar.pol)
		}
		if(mdlType == "ARMA22")
		{
			phis <- param[[3]][[i]][1:2]
			thetas <- param[[3]][[i]][3:4]
			ma.pol <- polymult(c(1,thetas),ma.pol)
			ar.pol <- polymult(c(1,-1*phis),ar.pol)
		}

		ma.null <- NULL
		if(length(ma.pol)==1) ma.null <- 0
		ma.pols[[length(ma.pols)+1]] <- c(ma.null,ma.pol)
		if(length(ma.pol)==1) ma.pols[[length(ma.pols)]] <- ma.pols[[length(ma.pols)]][-1]
		ar.null <- NULL
		if(length(ar.pol)==1) ar.null <- 0
		ar.pols[[length(ar.pols)+1]] <- c(ar.null,ar.pol)
		if(length(ar.pol)==1) ar.pols[[length(ar.pols)]] <- ar.pols[[length(ar.pols)]][-1]
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
			G.mat[i,j] <- ARMA2acf(ar = -1*ar.prod[-1],
				ma = ma.prod[-1],lag.max=1)[1]*ma.scale
		}
		g.acf <- ARMA2acf(ar = -1*ar.pols[[i]][-1],ma = ma.pols[[i]][-1],lag.max=Tdiff-1)
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

