mvar.lik <- function(x.acf,x)
{

	###########################
	#	mvar.lik
	#		by Tucker McElroy
	#
	#   Computes likelihood of a multivariate process
	#	via Levinson-Durbin algorithm
	#	x is the data, NxT matrix
	#	returns -2* log Gaussian likelihood, and resids
	#		
	##############################
	
	N <- dim(x)[1]
	T <- dim(x)[2]

	aseq <- solve(x.acf[,1,]) %*% x.acf[,2,]
	bseq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
	gamSeq <- NULL
	gamFlip <- NULL
	Qseq <- t(matrix(x[,1],ncol=1)) %*% solve(x.acf[,1,]) %*% 
		matrix(x[,1],ncol=1)
	logdet <- log(det(as.matrix(x.acf[,1,])))
	rhot <- matrix(0,nrow=N,ncol=N)
	Lam <- x.acf[,1,]
	Om <- x.acf[,1,]
	eps <- solve(chol(Lam)) %*% matrix(x[,1],ncol=1)
	for(t in 1:(T-2))
	{
		gamSeq <- cbind(x.acf[,t+1,],gamSeq)
		gamFlip <- rbind(gamFlip,x.acf[,t+1,])
		rhot <- gamSeq %*% aseq
		Lam <- x.acf[,1,] - gamSeq %*% bseq
		Om <- x.acf[,1,] - t(aseq) %*% gamFlip
		alphat <- t(bseq) %*% matrix(x[,1:t],ncol=1)
		new.eps <- solve(chol(Lam)) %*% (alphat - matrix(x[,t+1],ncol=1))
		eps <- rbind(eps,new.eps)
		Qseq <- Qseq + t(alphat - matrix(x[,t+1],ncol=1)) %*% 
			solve(Lam) %*% (alphat - matrix(x[,t+1],ncol=1))
		xit <- x.acf[,t+2,] - rhot
		bfact <- solve(Om) %*% t(xit)
		newb <- bseq - aseq %*% bfact
		afact <- solve(Lam) %*% xit
		newa <- aseq - bseq %*% afact
		bseq <- rbind(bfact,newb)
		aseq <- rbind(newa,afact)
		logdet <- logdet + log(det(Lam))
	} 
	gamSeq <- cbind(x.acf[,T,],gamSeq)
	Lam <- x.acf[,1,] - gamSeq %*% bseq
	logdet <- logdet + log(det(Lam))
	alphat <- t(bseq) %*% matrix(x[,1:(T-1)],ncol=1)
	new.eps <- solve(chol(Lam)) %*% (alphat - matrix(x[,T],ncol=1))
	eps <- rbind(eps,new.eps)
	eps <- matrix(eps,nrow=N,ncol=T)
	Qseq <- Qseq + t(alphat - matrix(x[,T],ncol=1)) %*% 
		solve(Lam) %*% (alphat - matrix(x[,T],ncol=1))
	lik <- Qseq + logdet 
	print(lik)

	return(list(c(Qseq,logdet),eps))

}
