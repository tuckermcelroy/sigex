specFactmvar <- function(xAcf)
{
	N <- dim(xAcf)[1]
	q <- dim(xAcf)[2] - 1
	Zzero <- matrix(0,nrow=N,ncol=N)
	eps <- 1
	thresh <- 10^(-15)
	Linvt <- diag(N)
	Dinv <- solve(xAcf[,1,])
	sqrtLam <- svd(xAcf[,1,])
	sqrtLam <- sqrtLam$u %*% diag(sqrt(sqrtLam$d),nrow=N) %*% t(sqrtLam$u)
	Gseq <- solve(sqrtLam)
	oldLam <- xAcf[,1,]
	Dsqrt <- sqrtLam
	gamSeq <- NULL
	count <- 1
	maxcount <- 500
	while((eps > thresh)&&(count<maxcount))
	{
		if(count <= q) nextgam <- xAcf[,count+1,] else nextgam <- Zzero
		gamSeq <- cbind(nextgam,gamSeq)
		Bseq <- gamSeq %*% Gseq
		Lam <- xAcf[,1,] - Bseq %*% t(Bseq)
		sqrtLam <- svd(Lam)
		sqrtLam <- sqrtLam$u %*% diag(sqrt(sqrtLam$d),nrow=N) %*% t(sqrtLam$u)
		Gseq <- rbind(cbind(Gseq,-1*Gseq%*%t(Bseq)%*%solve(sqrtLam)),
			cbind(t(rep(1,count) %x% Zzero),solve(sqrtLam)))
		count <- count+1
		if(count > q) eps <- sum((oldLam - Lam)^2)
		oldLam <- Lam
#		print(c(count,eps))
#		print(Lam)
	}
	Thetas <- cbind(Bseq,sqrtLam) %*% (diag(count) %x% solve(sqrtLam))
	arrayThetas <- array(Thetas,dim=c(N,N,count))
	out <- list(arrayThetas[,,(count-q):count],Lam,count)
	return(out)
}
