sigex.resid <- function(psi,mdl,data.ts,debug=FALSE)
{

	##########################################################################
	#
	#	sigex.resid
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
	#	Purpose: computes residuals from Gaussian likelihood of model
	#	Background:
	#		param is the name for the model parameters entered into
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Notes: handles missing values in data.ts, which are indicated by 1i
	#	Inputs:
	#		psi: see background.
	#		mdl: the specified sigex model, a list object
  #		data.ts: a T x N matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #   debug: set to TRUE if lik values should be printed to screen
  #	Outputs: a list with two elements
	#		resids: a N x (T-d) matrix of residuals, where d is
	#			the order of the full differencing polynomial.
	#			(A failure returns Inf.)
  #   x.acf: an array of dimension N x (T+1) x N of the autocovariances
  #     of the differenced model for lags 0 through T, corresponding to psi.
	#	Requires: sigex.zetalen, sigex.zeta2par, sigex.param2gcd, sigex.delta,
	#			mvar.midcast, sigex.acf
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	z <- x
	z[is.na(z)] <- 1i
	L.par <- mdl[[3]]
	D.par <- mdl[[3]]
	zeta.par <- vector("list",length(mdl[[3]]))
	acf.mat <- matrix(0,nrow=N*(T+1),ncol=N)

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
		acf.mat <- acf.mat + sigex.acf(L.par[[i]],D.par[[i]],mdl,i,zeta.par[[i]],delta,T+1)
	}

	x.acf <- array(acf.mat,dim=c(N,T+1,N))
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

	delta <- sigex.delta(mdl,0)
	attempt <- try(mvar.midcast(x.acf,z,delta,debug))
	if(!inherits(attempt, "try-error")) {
		lik.output <- attempt[[4]] } else lik.output <- Inf

	resids <- lik.output

	return(list(Re(resids),x.acf))
}


