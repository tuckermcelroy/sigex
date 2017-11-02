sigex.signal <- function(data,param,mdl,sigcomps)
{

	###############################
	#   sigex.signal
	#	by Tucker McElroy	
	#
	#	Computes signal extraction matrix and error cov matrix
	#		param must be in format yielded by sigex.default
	#		sigcomps provides indices of the desired components
	#     Output is in array form of dimension c(T,N,T,N), which are NxN number 
	#		of blocks each of size TxT, for filter matrix F and cov matrix V
	#
	#################################

	x <- t(data)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# compute differencing polynomials and matrices
	allcomps <- seq(1,length(mdl[[3]]))
	noisecomps <- allcomps[!allcomps %in% sigcomps]

	# get delta polynomials
	delta.data <- sigex.delta(mdl,0)
	delta.data <- array(t(delta.data %x% diag(N)),c(N,N,length(delta.data)))
	Tdiff <- T - dim(delta.data)[3] + 1
	delta.signal <- sigex.delta(mdl,noisecomps)
	delta.signal <- array(t(delta.signal %x% diag(N)),c(N,N,length(delta.signal)))	
	TdiffSig <- T - dim(delta.signal)[3] + 1
	delta.noise <- sigex.delta(mdl,sigcomps)
	delta.noise <- array(t(delta.noise %x% diag(N)),c(N,N,length(delta.noise)))	
	TdiffNoise <- T - dim(delta.noise)[3] + 1
#	print("HERE")

	delta.data2 <- array(0,c(N,N,T))
	delta.data2[,,1:dim(delta.data)[3]] <- delta.data[,,1:dim(delta.data)[3]]
	delta.data.mat <- sigex.blocktoep(delta.data2)[,dim(delta.data)[3]:T,,]
	delta.data.mat <- matrix(delta.data.mat,ncol=N*T)
#	print("HERE")

	delta.signal2 <- array(0,c(N,N,T))
	delta.signal2[,,1:dim(delta.signal)[3]] <- delta.signal[,,1:dim(delta.signal)[3]]
	delta.signal.bigmat <- sigex.blocktoep(delta.signal2)[,dim(delta.signal)[3]:T,,]
	delta.signal.bigmat <- matrix(delta.signal.bigmat,ncol=N*T)
	delta.signal.smallmat <- sigex.blocktoep(delta.signal2)[,dim(delta.data)[3]:T,,dim(delta.noise)[3]:T]
	delta.signal.smallmat <- matrix(delta.signal.smallmat,nrow=N*Tdiff)
#	print("HERE")

	delta.noise2 <- array(0,c(N,N,T))
	delta.noise2[,,1:dim(delta.noise)[3]] <- delta.noise[,,1:dim(delta.noise)[3]]
	delta.noise.bigmat <- sigex.blocktoep(delta.noise2)[,dim(delta.noise)[3]:T,,]
	delta.noise.bigmat <- matrix(delta.noise.bigmat,ncol=N*T)
	delta.noise.smallmat <- sigex.blocktoep(delta.noise2)[,dim(delta.data)[3]:T,,dim(delta.signal)[3]:T]
	delta.noise.smallmat <- matrix(delta.noise.smallmat,nrow=N*Tdiff)
#	print("HERE")

	######################################
	# compute autocovariances

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	acf.mat <- matrix(0,nrow=N*T,ncol=N)
	for(i in 1:length(mdl[[3]]))
	{
		L.par[[i]] <- param[[1]][[i]]
		D.par[[i]] <- param[[2]][[i]]
		mdlType <- mdl[[2]][i]	
		delta <- sigex.delta(mdl,i)
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,
			param[[3]][[i]],N,delta,T)		
	}
	x.acf <- array(acf.mat,dim=c(N,T,N))

	acfsignal.mat <- matrix(0,nrow=N*TdiffSig,ncol=N)
	for(i in sigcomps)
	{
		mdlType <- mdl[[2]][i]
		delta <- sigex.delta(mdl,c(noisecomps,i))
		acfsignal.mat <- acfsignal.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,
			param[[3]][[i]],N,delta,TdiffSig)		
	}
	signal.acf <- array(acfsignal.mat,dim=c(N,TdiffSig,N))

	acfnoise.mat <- matrix(0,nrow=N*TdiffNoise,ncol=N)
	for(i in noisecomps)
	{
		mdlType <- mdl[[2]][i]
		delta <- sigex.delta(mdl,c(sigcomps,i))
		acfnoise.mat <- acfnoise.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,
			param[[3]][[i]],N,delta,TdiffNoise)		
	}
	noise.acf <- array(acfnoise.mat,dim=c(N,TdiffNoise,N))

	######################################################
	# compute covariance matrices and inverses

	fulldiff <- sigex.delta(mdl,0)
	aseq <- solve(x.acf[,1,]) %*% x.acf[,2,]
	bseq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
	gamSeq <- NULL
	gamFlip <- NULL
	rhot <- matrix(0,nrow=N,ncol=N)
	Lam <- x.acf[,1,]
	Om <- x.acf[,1,]
	Gaminv <- solve(x.acf[,1,])
	for(t in 1:(Tdiff-2))
	{
		gamSeq <- cbind(x.acf[,t+1,],gamSeq)
		gamFlip <- rbind(gamFlip,x.acf[,t+1,])
		rhot <- gamSeq %*% aseq
		Lam <- x.acf[,1,] - (gamSeq %*% bseq + t(bseq) %*% t(gamSeq))/2
		Om <- x.acf[,1,] - (t(gamFlip) %*% aseq + t(aseq) %*% gamFlip)/2
		Ominv <- solve(Om)
		Gaminv <- rbind(cbind(Ominv, -1*Ominv %*% t(aseq)),
			cbind(-1*aseq %*% Ominv, Gaminv + aseq %*% Ominv %*% t(aseq)))
		xit <- x.acf[,t+2,] - rhot
		bfact <- solve(Om) %*% t(xit)
		newb <- bseq - aseq %*% bfact
		afact <- solve(Lam) %*% xit
		newa <- aseq - bseq %*% afact
		bseq <- rbind(bfact,newb)
		aseq <- rbind(newa,afact)
#		print(t)
	} 
	gamSeq <- cbind(x.acf[,Tdiff,],gamSeq)
	gamFlip <- rbind(gamFlip,x.acf[,Tdiff,])
	rhot <- gamSeq %*% aseq
	Lam <- x.acf[,1,] - (gamSeq %*% bseq + t(bseq) %*% t(gamSeq))/2
	Om <- x.acf[,1,] - (t(gamFlip) %*% aseq + t(aseq) %*% gamFlip)/2
	Ominv <- solve(Om)
	Gaminv <- rbind(cbind(Ominv, -1*Ominv %*% t(aseq)),
		cbind(-1*aseq %*% Ominv, Gaminv + aseq %*% Ominv %*% t(aseq)))

# Construct K matrix
#	Kcommut <- function(vect,m,n)
#	{
#		return(matrix(t(matrix(vect,nrow=m,ncol=n)),ncol=1))
#	}
#	temp <- apply(Gaminv,2,Kcommut,N,Tdiff)
#	KGaminv <- t(apply(t(temp),2,Kcommut,N,Tdiff))
#	GaminvArray <- array(KGaminv,dim=c(Tdiff,N,Tdiff,N))
	
	signal.acf2 <- array(0,c(dim(signal.acf)[1],dim(signal.acf)[3],dim(signal.acf)[2]))
	for(t in 1:dim(signal.acf2)[3]) { signal.acf2[,,t] <- signal.acf[,t,] }
	signal.acf2[,,1] <- signal.acf[,1,]/2
	signal.toep <- matrix(sigex.blocktoep(signal.acf2),N*(TdiffSig),N*(TdiffSig))
	signal.toep <- signal.toep + t(signal.toep)
#	print("HERE")

	noise.acf2 <- array(0,c(dim(noise.acf)[1],dim(noise.acf)[3],dim(noise.acf)[2]))
	for(t in 1:dim(noise.acf2)[3]) { noise.acf2[,,t] <- noise.acf[,t,] }
	noise.acf2[,,1] <- noise.acf[,1,]/2
	noise.toep <- matrix(sigex.blocktoep(noise.acf2),N*(TdiffNoise),N*(TdiffNoise))
	noise.toep <- noise.toep + t(noise.toep)
#	print("HERE")

	# compute sig ex formulas
	m.mat <- t(delta.signal.bigmat) %*% delta.signal.bigmat + t(delta.noise.bigmat) %*% delta.noise.bigmat
	p.mat <- t(delta.signal.bigmat) %*% signal.toep %*% t(delta.noise.smallmat) -
			t(delta.noise.bigmat) %*% noise.toep %*% t(delta.signal.smallmat)
	m.inv <- solve(m.mat)
	f.mat <- t(delta.noise.bigmat) %*% delta.noise.bigmat + p.mat %*% Gaminv %*% delta.data.mat
	f.mat <- m.inv %*% f.mat
	v.mat <- t(delta.noise.bigmat) %*% noise.toep %*% delta.noise.bigmat + 
			t(delta.signal.bigmat) %*% signal.toep %*% delta.signal.bigmat
	v.mat <- v.mat - p.mat %*% Gaminv %*% t(p.mat)
	v.mat <- m.inv %*% v.mat %*% m.inv

	return(list(f.mat,v.mat))

}

