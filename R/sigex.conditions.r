sigex.conditions <- function(data,psi,mdl)
{
	
	########################################
	#	SigEx conditions
	#		by Tucker McElroy
	#
	#	Computes condition number for each
	#		latent component innovation covariance matrix
	#		log determinant of correlation matrix
	#		(always -Inf for reduced rank models)
	#
	###########################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	psi <- Re(psi)

	# get xi portion
	ind <- 0
	A.mat <- matrix(0,N,N)
	A.mat[lower.tri(A.mat)] <- 1
	conds <- NULL
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
		k.order <- length(D.psi)
	
		cov.mat <- L.mat %*% diag(exp(D.psi),nrow=k.order) %*% t(L.mat)
		eta2 <- getGCD(cov.mat,N)[[2]]/diag(cov.mat)
		eta2[diag(cov.mat)==0] <- 0
		eta2[1] <- 1
		conds <- rbind(conds,eta2)	
	}

	return(conds)
}