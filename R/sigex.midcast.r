sigex.midcast <- function(psi,mdl,data.ts,castspan)
{

	##########################################################################
	#
	#	sigex.midcast
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
	#	Purpose: computes predictors for variables at various indices
	#	Background:	
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter manifold)
	#	Inputs:
	#		psi: see background. 
	#		mdl: the specified sigex model, a list object
  #		data.ts: a T x N matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #   castspan: an non-negative integer horizon giving number of fore- and aft-casts
  #	Outputs:
	#		list containing casts.x and casts.var 
	#		casts.x: N x H matrix of forecasts, midcasts, aftcasts, where H
	#			is the total number of time indices with missing values,
	#			given by cardinality( leads setminus {1,2,...,T} )	
	#		casts.var: NH x NH matrix of covariances of casting errors.
	#			note that casts.var.array <- array(casts.var,c(N,H,N,H)) 
	#			corresponds to cast.var.array[,j,,k] equal to the 
	#			covariance between the jth and kth casting errors
	#	Notes: presumes that regression effects have already been removed.
	#	Requires: sigex.param2gcd, sigex.zeta2par, sigex.zetalen, sigex.acf, sigex.delta,
	#			mvar.midcast, sigex.i2rag
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]
	psi <- Re(psi)

	z <- x
	z[is.na(z)] <- 1i
	out <- sigex.i2rag(z)
	leads.mid <- out[[1]]
  leads.fore <- NULL
  leads.aft <- NULL
	if(castspan > 0)
	{
	  leads.fore <- seq(T+1,T+castspan)
	  leads.aft <- seq(1-castspan,0)
	}
 
	#  alter z, inserting 1i for any out-of-sample forecasts/aftcasts
	if(length(leads.fore)>0) { z <- cbind(z,matrix(1i,nrow=N,ncol=length(leads.fore))) }
	if(length(leads.aft)>0) { z <- cbind(matrix(1i,nrow=N,ncol=length(leads.aft)),z) }
	TH <- dim(z)[2]

	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	zeta.par <- vector("list",length(mdl[[3]]))
	acf.mat <- matrix(0,nrow=N*TH,ncol=N)
	
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
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,zeta.par[[i]],delta,TH)		
	}

	x.acf <- array(acf.mat,dim=c(N,TH,N))
	reg.vec <- beta.par	

	# subtract regression effects from available sample only
	ind <- 0
	my.times <- seq(1+length(leads.aft),T+length(leads.aft))
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		z[k,my.times] <- z[k,my.times] - reg.mat %*% reg.vec[(ind+1):(ind+len)]
		ind <- ind+len
	}

	delta <- sigex.delta(mdl,0)
	attempt <- try(mvar.midcast(x.acf,z,delta,debug=FALSE))
	if(!inherits(attempt, "try-error")) {
		casts.x <- attempt[[1]]
		casts.var <- attempt[[2]] }
	
	return(list(casts.x,casts.var))
}


