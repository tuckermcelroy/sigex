sigex.whittle <- function(psi,mdl,data.ts)
{

	##########################################################################
	#
	#	sigex.whittle
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
	#	Purpose: computes Whittle likelihood of model 
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Notes: does not yet handle missing values in data.ts!!!
	#	Inputs:
	#		psi: see background. 
	#		mdl: the specified sigex model, a list object
	#		data.ts: a T x N matrix ts object; any missing values 
	#			must be encoded with NA in that entry
	#	Outputs:
	#		returns Whittle likelihood, evaluated at Fourier frequencies
	#			corresponding to model, for differenced time series,
	#			with hyper-parameter psi. 
	#	Requires: sigex.zetalen, sigex.zeta2par, sigex.param2gcd, sigex.delta,
	#			sigex.acf, sigex.spectra, sigex.psi2par
	#
	####################################################################
 
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]
	param <- sigex.psi2par(psi,mdl,data.ts)

	z <- x
	z[is.na(z)] <- 1i
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
		mdlType <- mdl[[2]][[i]]	
		delta <- mdl[[3]][[i]]
		zetalen <- sigex.zetalen(mdlType)
		if(zetalen > 0) {
			subzeta <- zeta[(ind+1):(ind+zetalen)]
			zeta.par[[i]] <- sigex.zeta2par(subzeta,mdlType)
		}
		ind <- ind + zetalen

		delta <- sigex.delta(mdl,i)
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,zeta.par[[i]],delta,T)		
	}

	x.acf <- array(acf.mat,dim=c(N,T,N))
	reg.vec <- beta.par	

	# subtract regression effects from available sample only
	ind <- 0
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		z[k,] <- z[k,] - reg.mat %*% reg.vec[(ind+1):(ind+len)]
		ind <- ind+len
	}

	# difference the data
	delta <- sigex.delta(mdl,0)
	x.diff <- as.matrix(filter(t(z),delta,method="convolution",
		sides=1)[length(delta):T,])
	Tdiff <- dim(x.diff)[1]
	x.diff <- t(x.diff)

	if(Tdiff %% 2 == 1) { grid <- Tdiff-1 } else { grid <- Tdiff-2 }
	T.m <- grid/2

	f.all <- t(rep(0,grid+1) %x% diag(N))
	for(i in 1:length(mdl[[3]]))
	{
		L.par <- param[[1]][[i]]
		D.par <- param[[2]][[i]]
		delta <- sigex.delta(mdl,i)
		f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],delta,grid)
		f.all <- f.all + matrix(f.comp,nrow=N)
	}
	f.all  <- array(f.all,c(N,N,(grid+1)))

	x.dft <- NULL
	for(i in 1:N)
	{
		x.dft <- rbind(x.dft,fft(x.diff[i,1:(grid+1)])*exp(-1i*(seq(1,grid+1)+T.m)/grid)/sqrt(grid+1))
	}
	
	whittle.lik <- 0
	for(l in 1:(grid+1))
	{
		whittle.lik <- whittle.lik + t(Conj(x.dft[,l])) %*% solve(f.all[,,l]) %*% x.dft[,l] +
			+ sum(log(Re(eigen(f.all[,,l])$values)))
	}
	whittle.lik <- Re(whittle.lik)/(grid+1)	
	print(whittle.lik)

	return(whittle.lik)

}
