mvar.cast <- function(x.acf,z)
{

	###########################
	#	mvar.cast
	#		by Tucker McElroy
	#
	#   Computes forecasts of a multivariate process
	#	via Levinson-Durbin algorithm
	#	z is the differenced (stationary) data, NxT matrix, with NA for various t,
	#	 which are to be imputed.
	#	presumes first observation is not all NA!
	#	NA or missing values, must be encoded with 1i in that entry
	#	returns casted tilde{x}, NxT, with imputations for the NA
	#		
	##############################
	
	N <- dim(z)[1]
	T <- dim(z)[2]
	all.series <- seq(1,N)
	all.indices <- seq(1,T)
	full.indices <- all.indices[colSums(z==1i)==0]
	cast.indices <- setdiff(all.indices,full.indices)

	cast.series <- all.series[z[,1]==1i]
	full.series <- setdiff(all.series,cast.series)
	full.mat <- diag(N)[full.series,]
	yhat <- x.acf[,1,] %*% t(full.mat) %*% solve(full.mat %*% x.acf[,1,] %*% t(full.mat)) %*%
		as.matrix(z[,1])
	y <- Re(yhat)
	aseq <- solve(x.acf[,1,]) %*% x.acf[,2,]
	bseq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
	gamSeq <- NULL
	gamFlip <- NULL
	rhot <- matrix(0,nrow=N,ncol=N)
	Lam <- x.acf[,1,]
	Om <- x.acf[,1,]
	for(t in 1:(T-2))
	{
		gamSeq <- cbind(x.acf[,t+1,],gamSeq)
		gamFlip <- rbind(gamFlip,x.acf[,t+1,])
		rhot <- gamSeq %*% aseq
		Lam <- x.acf[,1,] - gamSeq %*% bseq
		Om <- x.acf[,1,] - t(aseq) %*% gamFlip
		if((t+1) %in% cast.indices) 
		{
			cast.series <- all.series[z[,t+1]==1i]
			full.series <- setdiff(all.series,cast.series)
#			cast.mat <- diag(N)[cast.series,]
 			yhat <- t(bseq) %*% matrix(y[,1:t],ncol=1)
		 	if(length(full.series)>0) { 
				full.mat <- diag(N)[full.series,] 
				yhat <-  yhat +  Lam %*% t(full.mat) %*%
						solve(full.mat %*% Lam %*% t(full.mat)) %*% 
						(as.matrix(y[full.series,t]) - full.mat %*% yhat)	
			}
			y <- cbind(y,Re(yhat))
		} else { y <- cbind(y,Re(z[,t+1])) }
		xit <- x.acf[,t+2,] - rhot
		bfact <- solve(Om) %*% t(xit)
		newb <- bseq - aseq %*% bfact
		afact <- solve(Lam) %*% xit
		newa <- aseq - bseq %*% afact
		bseq <- rbind(bfact,newb)
		aseq <- rbind(newa,afact)
	} 
	gamSeq <- cbind(x.acf[,T,],gamSeq)
	Lam <- x.acf[,1,] - gamSeq %*% bseq
	if(T %in% cast.indices) 
	{
		cast.series <- all.series[z[,T]==1i]
		full.series <- setdiff(all.series,cast.series)
#		cast.mat <- diag(N)[cast.series,]
		yhat <- t(bseq) %*% matrix(y[,1:(T-1)],ncol=1)
	 	if(length(full.series)>0) { 
			full.mat <- diag(N)[full.series,] 
			yhat <- yhat + Lam %*% t(full.mat) %*%
					solve(full.mat %*% Lam %*% t(full.mat)) %*% 
					(as.matrix(y[full.series,T]) - full.mat %*% yhat)
		}
		y <- cbind(y,Re(yhat))
	} else { y <- cbind(y,Re(z[,T])) }

	return(y)
}
