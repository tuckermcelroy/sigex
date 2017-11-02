sigex.acf <- function(L.par,D.par,mdl,comp,mdlPar,N,delta,maxlag)
{
	
	###########################
	#	sigex.acf
	#		by Tucker McElroy
	#
	#   Computes acf of a latent multivariate component process
	#	L.par is the entries of unit lower triangular in GCD of
	#		the white noise covariance matrix
	#	D.par is the log entries of diagonal matrix in GCD
	#	mdl is the   model object
	#	mdlPar is that portion of param corresponding to mdlType
	#	N is the cross-sectional dimension
	#	comp is the component, and delta is a given
	#	    differencing polynomial; the acf of
	#		delta (B) x_t 
	#		is calculated at lags 0 through maxlag-1.
	#
	##############################

	mdlType <- mdl[[2]][comp]
	d.delta <- length(delta)
	xi.mat <- L.par %*% diag(exp(D.par),nrow=length(D.par)) %*% t(L.par)
	
	# get acf of stationary component
	if(mdlType == "wn") 
	{
		cycle.order <- 0
		rho <- 0
		omega <- 0 
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		cycle.MA <- out[[2]]
		psi.acf <- ARMAacvf(ar = -1*cycle.AR[-1],ma = polymul(cycle.MA,delta)[-1],
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
	}
	if(mdlType == "canonWN")
	{
		canon.delta <- mdl[[3]][[comp]]
		canon.d <- length(canon.delta)
		freq0 <- sum(canon.delta)^2 
		freqpi <- sum(canon.delta*((-1)^(seq(0,canon.d-1))))^2
		freqmax <- max(freq0,freqpi)								
		ma.poly <- polymul(rev(canon.delta),canon.delta)
		psi.acf <- c(1,rep(0,canon.d-1))
		psi.acf <- psi.acf - freqmax^{-1}*ma.poly[canon.d:(2*canon.d-1)]				
		psi.acf <- c(rev(psi.acf),psi.acf[-1])
		psi.ma <- Re(specFact(psi.acf))		
 		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymul(delta,psi.ma)
		psi.acf <- ARMAacvf(ar = NULL, ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]	
		psi.acf <- psi.acf*psi.scale
		x.acf <- psi.acf %x% xi.mat
	}
	if(mdlType == "AR1") 
	{
		cycle.order <- 1
		rho <- mdlPar[1]
		omega <- 0 
		out <- sigex.getcycle(cycle.order,rho,omega)
		cycle.AR <- out[[1]]
		cycle.MA <- out[[2]]
		psi.acf <- ARMAacvf(ar = -1*cycle.AR[-1],ma = polymul(cycle.MA,delta)[-1],
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
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
#		psi.acf <- ARMAacvf(ar = -1*cycle.AR[-1],ma = polymul(cycle.MA,delta)[-1],
#			lag.max=maxlag)[1:maxlag]
		psi.marep <- ARMAtoMA(ar = -1*cycle.AR[-1],ma = polymul(cycle.MA,delta)[-1],
			lag.max=maxlag)[1:maxlag]
		psi.acf <- ARMAacvf(ar = NULL,ma = psi.marep,lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
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
#		psi.acf <- ARMAacvf(ar = NULL,ma = cycle.MA[-1],lag.max=length(cycle.MA))
		psi.acf <- polymul(cycle.MA,rev(cycle.MA))
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
#		psi.acf <- c(rev(psi.acf),psi.acf[-1])
		psi.acf <- polysum(psi.acf,-freqmin*polymul(cycle.AR,rev(cycle.AR))) + 1e-10
	 	psi.ma <- Re(specFact(psi.acf))
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymul(delta,psi.ma)
		psi.AR <- cycle.AR
#		psi.acf <- ARMAacvf(ar = -1*psi.AR[-1],ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]
		psi.marep <- ARMAtoMA(ar = -1*psi.AR[-1],ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- ARMAacvf(ar = NULL,ma = psi.marep,lag.max=maxlag)[1:maxlag]
 		psi.acf <- psi.acf*psi.scale
	#	psi.acf <- filter(rev(filter(psi.acf,delta,method="convolution",sides=2)),
	#		delta,method="convolution",sides=2)
	#	psi.acf <- rev(psi.acf[d.delta:(maxlag + d.delta - 1)])
		x.acf <- psi.acf %x% xi.mat
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
		ma.acf <- c(rev(ma.acf),ma.acf[-1]) + 1e-10
		psi.ma <- Re(specFact(ma.acf))	
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymul(delta,psi.ma)
		psi.AR <- cycle.AR
#		psi.acf <- ARMAacvf(ar = -1*psi.AR[-1],ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]
		psi.marep <- ARMAtoMA(ar = -1*psi.AR[-1],ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- ARMAacvf(ar = NULL,ma = psi.marep,lag.max=maxlag)[1:maxlag]

#	Alternate method based on Trimbur (2010) formulas; our method uses ARMA factorization
#		alphas <- NULL
#		for(k in 1:cycle.order)
#		{
#			r <- seq(0,k-1)
#			alpha.sum <- sum(choose(k-1,r)*choose(cycle.order-1,r+cycle.order-k)*rho^(2*r))
#			new.alpha <- ((1-rho^2)^(cycle.order-k))*alpha.sum
#			alphas <- c(alphas,new.alpha)	
#		}
#		alphas <- c(rep(0,maxlag),alphas)
#		psi.acf <- rev(alphas)[1]
#		for(h in 1:(maxlag-1))
#		{
#			new.acf <- (rho^h)*cos(h*pi*omega)*sum(choose(h,seq(0,h))*rev(alphas)[1:(h+1)])
#			psi.acf <- c(psi.acf,new.acf)
#		}
#		psi.acf <- psi.acf*((1-rho^2)^(1-2*cycle.order))

		psi.acf <- psi.acf*psi.scale
		x.acf <- psi.acf %x% xi.mat
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
		ma.acf <- polysum(ma.acf,-freqmin*polymul(cycle.AR,rev(cycle.AR))) + 1e-10
		psi.ma <- Re(specFact(ma.acf))	
		psi.scale <- psi.ma[1]^2
		psi.ma <- psi.ma/psi.ma[1]	
		psi.MA <- polymul(delta,psi.ma)
		psi.AR <- cycle.AR
#		psi.acf <- ARMAacvf(ar = -1*psi.AR[-1],ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]
		psi.marep <- ARMAtoMA(ar = -1*psi.AR[-1],ma = psi.MA[-1],lag.max=maxlag)[1:maxlag]
		psi.acf <- ARMAacvf(ar = NULL,ma = psi.marep,lag.max=maxlag)[1:maxlag]
		psi.acf <- psi.acf*psi.scale
		x.acf <- psi.acf %x% xi.mat
	}
#	if(mdlType == "VAR1") 
#	{ 
#		phi.array <- array(mdlPar,c(N,N,1))
#		x.acf <- VARMAauto(phi.array,NULL,psi.mat,maxlag+d.delta-1)
#		x.acf <- matrix(x.acf,nrow=N)
#		x.acf <- t(matrix(apply(x.acf,3,t),nrow=N))	
#	}
	if(mdlType == "ARMA22")
	{
		ma.poly <- c(1,mdlPar[3:4])
		psi.acf <- ARMAacvf(ar = mdlPar[1:2], ma = polymul(ma.poly,delta)[-1], 
			lag.max=maxlag)[1:maxlag]
		x.acf <- psi.acf %x% xi.mat
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
		varma.theta <- array(t(c(delta,0)) %x% diag(N) + t(c(0,delta)) %x% vma.theta,c(N,N,d.delta+1))
		varma.acf <- VARMAauto(varma.phi,varma.theta,vma.sigma,T-1)
		x.array <- array(0,c(N,T,N))
      	for(i in 1:T) { x.array[,i,] <- varma.acf[,,i] }
		x.acf <- matrix(x.array,c(N*T,T))
	}
	
	return(x.acf)
}
