sigex.par2psi <- function(param,flag,mdl)
{

	#################################
	#   sigex.par2psi
	#	by Tucker McElroy	
	#
	#	Takes list parameter param, with flag for fixed values, and produces
	#		full parameter vector (complex) psi as the inverse
	#		of sigex.psi2par
	#		psi = [xi,zeta,beta]
	#			xi ~ all pre-parameters for covariance matrices
	#			zeta ~ all pre-parameters for t.s. models
	#			beta ~ all regression parameters
	#	Also produces the flag vector:
	#		1 if the corresponding component is real
	#		0 if the corresponding component has an imaginary part
	#
	#################################

	xi <- NULL
	zeta <- NULL
	boundlist <- mdl[[5]]
	for(i in 1:length(mdl[[3]]))
	{
		bounds <- boundlist[[i]]
		mdlType <- mdl[[2]][i]
		delta <- mdl[[3]][[i]]
		vrank <- mdl[[1]][[i]]
		L.mat <- param[[1]][[i]]
		D.mat <- param[[2]][[i]]
		new.xi <- c(L.mat[lower.tri(diag(N))[,as.vector(vrank)]],D.mat)
		xi <- c(xi,new.xi)
		new.zeta <- sigex.par2zeta(param[[3]][[i]],mdlType,delta,bounds)
		zeta <- c(zeta,new.zeta)
	}

	beta <- param[[4]]
	psi <- c(xi,zeta,beta)
	psi <- psi + 1i*flag

	return(psi)
}
