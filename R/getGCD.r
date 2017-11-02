getGCD <- function(Sigma,Rank)
{

	#############################
	#	getGCD
	#		by Tucker McElroy
	#	Gets Generalized Cholesky Decomposition of Sigma,
	#	 allowing for zero Schur complements
	#      Rank is the presumed rank of the matrix, less than or equal
	#	  to the dimension of the matrix
	#	Output consists of the lower cholesky matrix,
	#	  and the diagonal matrix, of reduced dimension
	#
	#############################

	N <- dim(Sigma)[1]
	L.mat <- matrix(1,1,1)
	L.mat.inv <- L.mat
	D.mat <- Sigma[1,1]
	if(N > 1) {
	for(j in 2:N)
	{
		
		D.inv <- 1/D.mat
		D.inv[D.mat==0] <- 0
		new.sigma <- Sigma[j,1:(j-1)]
		if(j==2) { new.sigma <- as.matrix(new.sigma); L.mat <- as.matrix(L.mat) }
		new.l <- new.sigma %*% t(L.mat.inv)*D.inv
		new.l.tilde <- new.l %*% L.mat.inv
		L.mat <- cbind(L.mat,rep(0,(j-1)))
		L.mat <- rbind(L.mat,c(new.l,1))
		L.mat.inv <- cbind(L.mat.inv,rep(0,j-1))
		L.mat.inv <- rbind(L.mat.inv,c(-1*new.l.tilde,1))
		if(j==2) new.d <- Sigma[2,2] - new.l^2*D.mat
		if(j > 2) new.d <- Sigma[j,j] - new.l %*% diag(D.mat) %*% t(new.l)
		if(new.d <= 0) { new.d <- 0 }
		D.mat <- c(D.mat,new.d)
	} }
	
	rank.index <- rank(D.mat,ties.method="first")
	dims <- seq(1,N)[rank.index > (N-Rank)]

	L.mat <- matrix(L.mat[,dims],nrow=N,ncol=length(dims))
	D.mat <- D.mat[dims]

#	print(Lmat)
#	print(Dmat)
#	print(Lmat %*% diag(Dmat,nrow=length(dims)) %*% t(Lmat))

	return(list(L.mat,D.mat))
#	return(list(L.mat,D.mat,dims))
}



