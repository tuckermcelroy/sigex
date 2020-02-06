sigex.sim <- function(psi,mdl,simlen,burnin,dof,init)
{

	##########################################################################
	#
	#	sigex.sim
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

	################# Documentation #####################################
	#
	#	Purpose: simulate a stochastic process
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		psi: see background. 
	#		mdl: the specified sigex model, a list object
	#		simlen: length of the simulation
	#		burnin: initial stretch of simulation, later discarded
	#		dof: innovations of simulation are student t with dof
	#			degrees of freedom, set dof=Inf to get Gaussian
	#		init: initial values for process, should have length d,
	#			where d is the order of the full differencing polynomial
	#	Outputs:
	#		sims: matrix of dimension simlen x N of simulated stochastic process
	#			corresponding to mdl with parameter psi.
	#	Requires: sigex.delta, sigex.param2gcd, sigex.zeta2par, sigex.zetalen,
	#			sigex.acf
	#
	####################################################################

	N <- length(mdl[[4]])
	boundlist <- mdl[[5]]
	delta <- sigex.delta(mdl,0)
	d <- length(delta) - 1
	T <- simlen + burnin + d

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
			zeta.par[[i]] <- sigex.zeta2par(subzeta,mdlType,N,bounds)
		}
		ind <- ind + zetalen
	
		delta <- sigex.delta(mdl,i)
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,zeta.par[[i]],delta,T)		
	}

	x.acf <- array(acf.mat,dim=c(N,T,N))
	reg.vec <- beta.par	

	if(dof == Inf) { eps <- t(matrix(rnorm(N*T),nrow=N)) } else {
		eps <- t(matrix(rt(N*T,df=dof),nrow=N)) }
	aseq <- solve(x.acf[,1,]) %*% x.acf[,2,]
	bseq <- solve(x.acf[,1,]) %*% t(x.acf[,2,])
	gamSeq <- NULL
	gamFlip <- NULL
	rhot <- matrix(0,nrow=N,ncol=N)
	Lam <- x.acf[,1,]
	Om <- x.acf[,1,]
	sim <- chol(Lam) %*% matrix(eps[1,],ncol=1)
	for(t in 1:(T-2))
	{
		gamSeq <- cbind(x.acf[,t+1,],gamSeq)
		gamFlip <- rbind(gamFlip,x.acf[,t+1,])
		rhot <- gamSeq %*% aseq
		Lam <- x.acf[,1,] - gamSeq %*% bseq
		Om <- x.acf[,1,] - t(aseq) %*% gamFlip
		alphat <- t(bseq) %*% matrix(sim,ncol=1)
		new.sim <- chol(Lam) %*% matrix(eps[(t+1),],ncol=1) + alphat
		sim <- rbind(sim,new.sim)
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
	alphat <- t(bseq) %*% matrix(sim,ncol=1)
	new.sim <- chol(Lam) %*% matrix(eps[T,],ncol=1) + alphat
	sim <- rbind(sim,new.sim)
	sim <- matrix(sim,nrow=N)
	delta <- sigex.delta(mdl,0)
	delta.recurse <- -delta[-1]/delta[1]
	
	sims <- as.matrix(filter(t(sim),delta.recurse,method="recursive",init)[-seq(1,d),])
	sims <- as.matrix(sims[(burnin+1):(burnin+simlen),])

	# add fixed effects here, regressors must have length simlen

	return(sims)
}


