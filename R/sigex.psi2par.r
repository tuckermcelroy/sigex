sigex.psi2par <- function(psi,mdl,data)
{

	#################################
	#   sigex.psi2par
	#	by Tucker McElroy	
	#
	#	Takes full parameter vector psi, and produces param format,
	#		corresponding L and D matrices,
	#		as well as model and regression parameters
	#		psi = [xi,zeta,beta]
	#			xi ~ all pre-parameters for covariance matrices
	#			zeta ~ all pre-parameters for t.s. models
	#			beta ~ all regression parameters
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	zeta.par <- vector("list",length(mdl[[3]]))
	psi <- Re(psi)
	boundlist <- mdl[[5]]

	# get xi portion
	ind <- 0
	A.mat <- matrix(0,N,N)
	A.mat[lower.tri(A.mat)] <- 1
	for(i in 1:length(mdl[[3]]))
	{
		vrank <- mdl[[1]][[i]]
		D.dim <- length(vrank)
		L.dim <- sum(A.mat[,as.vector(vrank)])
		L.psi <- NULL
		if(L.dim > 0) L.psi <- psi[(ind+1):(ind+L.dim)]
		ind <- ind+L.dim
		D.psi <- psi[(ind+1):(ind+D.dim)]
		ind <- ind+D.dim
		L.mat <- sigex.param2gcd(L.psi,N,as.vector(vrank))
		L.par[[i]] <- L.mat
		D.par[[i]] <- D.psi
	}
	
	# get beta portion
	beta.len <- 0
	for(i in 1:N) 
	{
		beta.len <- beta.len + dim(mdl[[4]][[i]])[2]
	}
	beta.par <- as.vector(psi[(length(psi)-beta.len+1):length(psi)])

	# get zeta portion
	if(length(psi)-beta.len-ind > 0) {
		zeta <- psi[(ind+1):(length(psi)-beta.len)] }
	ind <- 0
	for(i in 1:length(mdl[[3]]))
	{
		bounds <- boundlist[[i]]
		mdlType <- mdl[[2]][i]	
		delta <- mdl[[3]][[i]]
		zetalen <- sigex.zetalen(mdlType)
		if(zetalen > 0) {
			subzeta <- zeta[(ind+1):(ind+zetalen)]
			zeta.par[[i]] <- sigex.zeta2par(subzeta,mdlType,delta,N,bounds)
		}
		ind <- ind + zetalen
	} 

	return(list(L.par,D.par,zeta.par,beta.par))
}

