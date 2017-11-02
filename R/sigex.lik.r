sigex.lik <- function(psi,mdl,data)
{

	###########################
	#	sigex.lik
	#		by Tucker McElroy
	#
	#   Computes likelihood of a latent multivariate component process
	#		psi = [xi,zeta,beta]
	#			xi ~ all pre-parameters for covariance matrices
	#			zeta ~ all pre-parameters for t.s. models
	#			beta ~ all regression parameters
	#	mdl is the model structure
	#	x is the data, NxT matrix
	#	returns -2* log Gaussian likelihood
	#		
	##############################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]
	psi <- Re(psi)
	boundlist <- mdl[[5]]
 
	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	zeta.par <- vector("list",length(mdl[[3]]))
	acf.mat <- matrix(0,nrow=N*T,ncol=N)
	
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

		delta <- sigex.delta(mdl,i)
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,
			zeta.par[[i]],N,delta,T)		
	}

	x.acf <- array(acf.mat,dim=c(N,T,N))
	reg.vec <- beta.par	

	# subtract regression effects
	ind <- 0
	data.diff <- data
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		data.diff[,k] <- data.diff[,k] - reg.mat %*% reg.vec[(ind+1):(ind+len)]
		ind <- ind+len
	}

	# difference the data
	fulldiff <-  sigex.delta(mdl,0)
	x.diff <- as.matrix(filter(data.diff,fulldiff,method="convolution",
		sides=1))[length(fulldiff):T,] 
	Tdiff <- dim(x.diff)[1]
	x.diff <- t(x.diff)

	attempt <- try(mvar.lik(x.acf,x.diff),TRUE)
	if(!inherits(attempt, "try-error")) {
		lik.output <- attempt[[1]] } else lik.output <- Inf

#	lik.output <- mvar.lik(x.acf,x.diff)[[1]]
	
	return(sum(lik.output))
}
