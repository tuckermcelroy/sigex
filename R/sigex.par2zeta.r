sigex.par2zeta <- function(mdlPar,mdlType,delta,bounds)
{

	#################################
	#   sigex.par2zeta
	#	by Tucker McElroy	
	#
	#	Takes mdlPar as param portion for mdlType, and produces
	#		full parameter vector zeta 
	#		psi = [xi,zeta,beta]
	#			xi ~ all pre-parameters for covariance matrices
	#			zeta ~ all pre-parameters for t.s. models
	#			beta ~ all regression parameters
	#	Also produces the flag vector:
	#		1 if the corresponding component is real
	#		0 if the corresponding component has an imaginary part
	#	bounds gives bounds for rho and omega:
	#			rho in (bounds[1],bounds[2])
	#			omega in (bounds[3],bounds[4])
	#
	#################################

	# defaults
	low.rho <- bounds[1]
	upp.rho <- bounds[2]
	low.omega <- bounds[3]
	upp.omega <- bounds[4]
		
	if(mdlType %in% c("wn","canonWN")) { zeta <- NULL; zeta.flag <- NULL } 
	if(mdlType == "AR1") 
	{ 
		phi <- mdlPar[1]
		zeta <- phi
		zeta <- log(1 + zeta) - log(1 - zeta)
	}
	if(mdlType %in% c("cycleBW1","cycleBW2","cycleBW3","cycleBW4","cycleBW5",
		"cycleBW6","cycleBW7","cycleBW8","cycleBW9","cycleBW10",
		"canonCycleBW1","canonCycleBW2","canonCycleBW3","canonCycleBW4",
		"canonCycleBW5","canonCycleBW6","canonCycleBW7","canonCycleBW8",
		"canonCycleBW9","canonCycleBW10","cycleBAL1","cycleBAL2","cycleBAL3",
		"cycleBAL4","cycleBAL5","cycleBAL6","cycleBAL7","cycleBAL8",
		"cycleBAL9","cycleBAL10","canonCycleBAL1","canonCycleBAL2",
		"canonCycleBAL3","canonCycleBAL4","canonCycleBAL5","canonCycleBAL6",
		"canonCycleBAL7","canonCycleBAL8","canonCycleBAL9","canonCycleBAL10"))
	{
		x.rho <- (mdlPar[1] - low.rho)/(upp.rho - low.rho)
		rho <- log(x.rho) - log(1 - x.rho)
		x.omega <- (mdlPar[2] - low.omega)/(upp.omega - low.omega)
		omega <- log(x.omega) - log(1 - x.omega)
		zeta <- c(rho,omega)
	}
#	if(mdlType == "VAR1") 
#	{ 
#		phi1.mat <- param[[3]][[i]]
#		subphi <- Re(phi1.mat)
#		subphi <- as.vector(matrix(subphi,ncol=1))
#		subflag <- Im(phi1.mat)
#	}
	if(mdlType == "ARMA22")
	{
		phis <- mdlPar[1:2]
		thetas <- -1*mdlPar[3:4]
		rhos <- c(phis[1]/(1-phis[2]),phis[1]^2/(1-phis[2]) + phis[2])
		rhos[2] <- log(rhos[2]-2*rhos[1]^2+1) - log(1-rhos[2])
		rhos[1] <- log(1+rhos[1]) - log(1-rhos[1])
		zeta <- rhos
		rhos <- c(thetas[1]/(1-thetas[2]),thetas[1]^2/(1-thetas[2]) + thetas[2])
		rhos[2] <- log(rhos[2]-2*rhos[1]^2+1) - log(1-rhos[2])
		rhos[1] <- log(1+rhos[1]) - log(1-rhos[1])
		zeta <- c(zeta,rhos)
	}
	if(mdlType == "damped")
	{
		L.zeta <- mdlPar[[1]]
		D.zeta <- mdlPar[[2]]
		phi <- mdlPar[[3]]
		N <- dim(L.zeta)[1]
		L.zeta <- L.zeta[lower.tri(diag(N))]
		zeta <- c(L.zeta,D.zeta)
		phi <- log(1 + phi) - log(1 - phi)
		zeta <- c(zeta,phi)
	}

	return(zeta)
}
