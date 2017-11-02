sigex.default <- function(mdl,data)
{

	###############################
	#   sigex.default
	#	by Tucker McElroy	
	#
	#	Determines default parameter setting of all zeroes,
	#		for the given specified mdl
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	A.mat <- matrix(0,N,N)
	A.mat[lower.tri(A.mat)] <- 1
	psi.len <- 0
	boundlist <- mdl[[5]]

	for(i in 1:length(mdl[[3]]))
	{	
		vrank <- mdl[[1]][[i]]
		D.dim <- length(vrank)
		L.dim <- sum(A.mat[,as.vector(vrank)])
		psi.len <- psi.len + D.dim + L.dim
		mdlType <- mdl[[2]][i]
		psi.len <- psi.len + sigex.zetalen(mdlType)
	}
	for(k in 1:N)
	{
		psi.len <- psi.len + dim(mdl[[4]][[k]])[2]
	}
	psi <- rep(0,psi.len) + 1i*rep(1,psi.len)
	par.default <- sigex.psi2par(psi,mdl,data)
	flag <- Im(psi)

	return(list(par.default,flag))
}

