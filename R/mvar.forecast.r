mvar.forecast <- function(x.acf,z)
{

	###########################
	#	mvar.forecast
	#		by Tucker McElroy
	#
	#   Computes multi-step forecasts and predictors of a multivariate process
	#	via Levinson-Durbin algorithm
	#	z is the differenced (stationary) data, Nx(T+H) matrix, with NA for various t,
	#	 which are to be imputed.
	#	presumes first T observations are not NA, and latter H observations are NA
	#	NA or missing values, must be encoded with 1i in that entry
	#	returns casted tilde{x}, Nx(T+H), with imputations for the NA
	#	 also returns NTxH matrix of predictors
	#		
	##############################
	
	N <- dim(z)[1]
	TH <- dim(z)[2]
	all.series <- seq(1,N)
	all.indices <- seq(1,TH)
	full.indices <- all.indices[colSums(z==1i)==0]
	cast.indices <- setdiff(all.indices,full.indices)
	H <- length(cast.indices)
	T <- length(full.indices)

#	pred.stack <- diag(N*T)
	pred.stack <- NULL
	aseq <- solve(x.acf[,1,]) %*% x.acf[,2,]
	bseq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
	gamSeq <- NULL
	gamFlip <- NULL
	rhot <- matrix(0,nrow=N,ncol=N)
	Lam <- x.acf[,1,]
	Om <- x.acf[,1,]
	for(t in 1:(TH-2))
	{
		gamSeq <- cbind(x.acf[,t+1,],gamSeq)
		gamFlip <- rbind(gamFlip,x.acf[,t+1,])
		rhot <- gamSeq %*% aseq
		Lam <- x.acf[,1,] - gamSeq %*% bseq
		Om <- x.acf[,1,] - t(aseq) %*% gamFlip
		xit <- x.acf[,t+2,] - rhot
		bfact <- solve(Om) %*% t(xit)
		newb <- bseq - aseq %*% bfact
		afact <- solve(Lam) %*% xit
		newa <- aseq - bseq %*% afact
		bseq <- rbind(bfact,newb)
		aseq <- rbind(newa,afact)
		if(t==(T-1))
		{
			pred.next <- bseq
			pred.stack <- pred.next
		}
		if(t > (T-1)) 
		{
			bseq.fore <- as.matrix(bseq[1:(N*T),])
			bseq.aft <- as.matrix(bseq[(N*T +1):((t+1)*N),]) 
			pred.next <- bseq.fore + pred.stack %*% bseq.aft 		
			pred.stack <- cbind(pred.stack,pred.next)
		}
	}
	x.cast <- t(pred.stack) %*% matrix(Re(z[,1:T]),ncol=1)
	x.cast <- matrix(x.cast,nrow=N)
	y.full <- cbind(z[,1:T],x.cast)
	
	return(list(y.full,pred.stack)) 
}
