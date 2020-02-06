sigex.conditions <- function(data.ts,psi,mdl)
{
	
	##########################################################################
	#
	#	sigex.conditions
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
	#	Purpose: computes condition number for a covariance matrix
	#	Background: a non-negative definite matrix Sigma has a
	#		Generalized Cholesky Decomposition (GCD) of the form
	#		Sigma = L %*% D %*% t(L),
	#		where L is unit lower triangular and D is diagonal with
	#		non-negative entries, referred to as the Schur complements
	#		of Sigma.  The number of nonzero Schur complements equals
	#		the rank of Sigma.  The condition numbers can be computed
	#		by dividing D by the diagonal of Sigma.
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		psi: see background.  
	#		mdl: the specified sigex model, a list object
	#	Outputs:
	#		conds: a S x N matrix of condition numbers, where
	#			S is the number of components.  Each row gives
	#			the N condition numbers for the innovation
	#			covariance matrix of the corresponding latent component.
	#	Requires: sigex.param2gcd, getGCD
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# get xi portion
	ind <- 0
	A.mat <- matrix(0,N,N)
	A.mat[lower.tri(A.mat)] <- 1
	conds <- NULL
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
		k.order <- length(D.psi)
	
		cov.mat <- L.mat %*% diag(exp(D.psi),nrow=k.order) %*% t(L.mat)
		eta2 <- getGCD(cov.mat,N)[[2]]/diag(cov.mat)
		eta2[diag(cov.mat)==0] <- 0
		eta2[1] <- 1
		conds <- rbind(conds,eta2)	
	}

	return(conds)
}