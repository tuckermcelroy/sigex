sigex.cast <- function(psi,mdl,data.ts,leads)
{

	##########################################################################
	#
	#	sigex.cast
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
	#	Purpose: computes forecasts and aftcasts, without uncertainty;
	#		a faster version of sigex.midcast, if you don't need midcasts
	#	Background:	
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		psi: see background. 
	#		mdl: the specified sigex model, a list object
	#		data.ts: a T x N matrix ts object (with no missing values)
	#			corresponding to N time series of length T
	#		leads: an integer sequence of desired casts; include index t
	#			to obtain an estimate of x_t.  These integers don't have
	#			to be a subset of {1,2,...,T}.  Include integers greater than
	#			T to get forecasts, or less than 1 to get aftcasts. 
	#	Outputs:
	#		x.casted: N x (T+H) matrix of aftcasts, data, and forecasts, where H
	#			is the total number of forecasts and aftcasts	
	#	Notes: presumes that regression effects have already been removed.
	#	Requires: sigex.param2gcd, sigex.zeta2par, sigex.zetalen, sigex.acf, sigex.delta,
	#			mvar.forecast
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	indices <- union(seq(1,T),leads)
	aft.index <- min(indices)
	fore.index <- max(indices)

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	zeta.par <- vector("list",length(mdl[[3]]))
	acf.mat <- matrix(0,nrow=N*length(indices),ncol=N)
	
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
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,zeta.par[[i]],delta,length(indices))		
	}

	x.acf <- array(acf.mat,dim=c(N,length(indices),N))
	reg.vec <- beta.par	

	# subtract regression effects
	ind <- 0
	data.diff <- data.ts
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		data.diff[,k] <- data.diff[,k] - reg.mat %*% reg.vec[(ind+1):(ind+len)]
		ind <- ind+len
	}

	# difference the data
	fulldiff <-  sigex.delta(mdl,0)
	del <- length(fulldiff) - 1
	x.diff <- as.matrix(filter(data.diff,fulldiff,method="convolution",
		sides=1)[length(fulldiff):T,])
	Tdiff <- dim(x.diff)[1]
	x.diff <- t(x.diff)
 
	fore.cast <- NULL
	aft.cast <- NULL
	if(fore.index > T)
	{
		x.fore <- cbind(x.diff,matrix(1i,N,(fore.index-T)))
		diff.cast <- mvar.forecast(x.acf,x.fore,FALSE)[[1]]
		if(del > 0) {
			fore.cast <- as.matrix(filter(init = matrix(data.diff[del:1,],ncol=N),
				x=t(diff.cast)/fulldiff[1],filter=-1*fulldiff[-1]/fulldiff[1],
				method="recursive"))
		} else { fore.cast <- t(diff.cast) }
		fore.cast <- as.matrix(fore.cast[(Tdiff+1):(Tdiff+fore.index-T),])
		fore.cast <- t(fore.cast)		
	}
	if(aft.index < 1)
	{
		x.rev <- t(as.matrix(t(x.diff)[seq(Tdiff,1),]))
		x.aft <- cbind(x.rev,matrix(1i,N,(1-aft.index)))
		diff.cast <- mvar.forecast(aperm(x.acf,c(3,2,1)),x.aft,FALSE)[[1]]
 		if(del > 0) {
			aft.cast <- as.matrix(filter(init = matrix(data.diff[(Tdiff+1):T,],ncol=N),
				x=t(diff.cast)/fulldiff[del+1],filter=-1*rev(fulldiff)[-1]/fulldiff[del+1],
				method="recursive"))
		} else { aft.cast <- t(diff.cast) }
		aft.cast <- as.matrix(aft.cast[(Tdiff+1):(Tdiff+1-aft.index),])
		aft.cast <- as.matrix(aft.cast[seq(1-aft.index,1),])
		aft.cast <- t(aft.cast)		
	}
	x.casted <- cbind(aft.cast,t(data.diff),fore.cast)
	x.real <- x.casted
	if(length(aft.cast) > 0) { x.real[,1:dim(aft.cast)[2]] <- NA }
	if(length(fore.cast) > 0) { x.real[,(dim(x.real)[2]-dim(fore.cast)[2]):dim(x.real)[2]] <- NA }

	return(x.casted)
}


