sigex.signal <- function(data.ts,param,mdl,sigcomps)
{

	##########################################################################
	#
	#	sigex.signal
	# 	    Copyright (C) 2017  Tucker McElroy
	#
	#    This program is free software: you can redistribute it and/or modify
	#    it under the terms of the GNU General Public License as published by
	#    the Free Software Foundation, either version 3 of the License, or
	#    (at your option) any later version.
	#
	#    This program is distributed in the hope that it will be useful,
	#    but WITHOUT ANY WARRANTY; without even the implied warranty of
	#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	#    GNU General Public License for more details.
	#
	#    You should have received a copy of the GNU General Public License
	#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
	#
	############################################################################

	################# Documentation ############################################
	#
	#	Purpose: computes signal extraction matrix and error covariance matrix
	#	Background:	
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		We have a model for each w_t process, and can compute its
	#		autocovariance function (acf), and denote its autocovariance
	#		generating function (acgf) via gamma_w (B).
	#		The signal extraction filter for y_t is determined from
	#		this acgf and delta.
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background.  Must have form specified by mdl
	#		mdl: the specified sigex model, a list object
	#		sigcomps: indices of the latent components composing the signal
	#	Outputs:
	#		list object of f.mat and v.mat
	#		f.mat: array of dimension c(T,N,T,N), where f.mat[,j,,k]
	#			is the signal extraction matrix that utilizes input series k
	#			to generate the signal estimate for series j.
	#		v.mat:  array of dimension c(T,N,T,N), where v.mat[,j,,k]
	#			is the error covariance matrix arising from input series k
	#			used to generate the signal estimate for series j.
	#	Requires: sigex.delta, sigex.blocktoep, sigex.acf
	#
	####################################################################
 
	x <- t(data.ts)
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

	delta.data2 <- array(0,c(N,N,T))
	delta.data2[,,1:dim(delta.data)[3]] <- delta.data[,,1:dim(delta.data)[3]]
	delta.data.mat <- sigex.blocktoep(delta.data2)[,dim(delta.data)[3]:T,,]
	delta.data.mat <- matrix(delta.data.mat,ncol=N*T)
 
	delta.signal2 <- array(0,c(N,N,T))
	delta.signal2[,,1:dim(delta.signal)[3]] <- delta.signal[,,1:dim(delta.signal)[3]]
	delta.signal.bigmat <- sigex.blocktoep(delta.signal2)[,dim(delta.signal)[3]:T,,]
	delta.signal.bigmat <- matrix(delta.signal.bigmat,ncol=N*T)
	delta.signal.smallmat <- sigex.blocktoep(delta.signal2)[,dim(delta.data)[3]:T,,dim(delta.noise)[3]:T]
	delta.signal.smallmat <- matrix(delta.signal.smallmat,nrow=N*Tdiff)
 
	delta.noise2 <- array(0,c(N,N,T))
	delta.noise2[,,1:dim(delta.noise)[3]] <- delta.noise[,,1:dim(delta.noise)[3]]
	delta.noise.bigmat <- sigex.blocktoep(delta.noise2)[,dim(delta.noise)[3]:T,,]
	delta.noise.bigmat <- matrix(delta.noise.bigmat,ncol=N*T)
	delta.noise.smallmat <- sigex.blocktoep(delta.noise2)[,dim(delta.data)[3]:T,,dim(delta.signal)[3]:T]
	delta.noise.smallmat <- matrix(delta.noise.smallmat,nrow=N*Tdiff)
 
	######################################
	# compute autocovariances

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	acf.mat <- matrix(0,nrow=N*T,ncol=N)
	for(i in 1:length(mdl[[3]]))
	{
		L.par[[i]] <- param[[1]][[i]]
		D.par[[i]] <- param[[2]][[i]]
		delta <- sigex.delta(mdl,i)
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,param[[3]][[i]],delta,T)		
	}
	x.acf <- array(acf.mat,dim=c(N,T,N))

	acfsignal.mat <- matrix(0,nrow=N*TdiffSig,ncol=N)
	for(i in sigcomps)
	{
		delta <- sigex.delta(mdl,c(noisecomps,i))
		acfsignal.mat <- acfsignal.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,param[[3]][[i]],delta,TdiffSig)		
	}
	signal.acf <- array(acfsignal.mat,dim=c(N,TdiffSig,N))

	acfnoise.mat <- matrix(0,nrow=N*TdiffNoise,ncol=N)
	for(i in noisecomps)
	{
		delta <- sigex.delta(mdl,c(sigcomps,i))
		acfnoise.mat <- acfnoise.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,param[[3]][[i]],delta,TdiffNoise)		
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
	} 
	gamSeq <- cbind(x.acf[,Tdiff,],gamSeq)
	gamFlip <- rbind(gamFlip,x.acf[,Tdiff,])
	rhot <- gamSeq %*% aseq
	Lam <- x.acf[,1,] - (gamSeq %*% bseq + t(bseq) %*% t(gamSeq))/2
	Om <- x.acf[,1,] - (t(gamFlip) %*% aseq + t(aseq) %*% gamFlip)/2
	Ominv <- solve(Om)
	Gaminv <- rbind(cbind(Ominv, -1*Ominv %*% t(aseq)),
		cbind(-1*aseq %*% Ominv, Gaminv + aseq %*% Ominv %*% t(aseq)))
	
	signal.acf2 <- array(0,c(dim(signal.acf)[1],dim(signal.acf)[3],dim(signal.acf)[2]))
	signal.acf2 <- aperm(signal.acf,c(1,3,2))
	signal.acf2[,,1] <- signal.acf[,1,]/2
	signal.toep <- matrix(sigex.blocktoep(signal.acf2),N*(TdiffSig),N*(TdiffSig))
	signal.toep <- signal.toep + t(signal.toep)

	noise.acf2 <- array(0,c(dim(noise.acf)[1],dim(noise.acf)[3],dim(noise.acf)[2]))
	noise.acf2 <- aperm(noise.acf,c(1,3,2))
	noise.acf2[,,1] <- noise.acf[,1,]/2
	noise.toep <- matrix(sigex.blocktoep(noise.acf2),N*(TdiffNoise),N*(TdiffNoise))
	noise.toep <- noise.toep + t(noise.toep)

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

