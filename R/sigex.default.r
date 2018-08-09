sigex.default <- function(mdl,data.ts)
{

	##########################################################################
	#
	#	sigex.default
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
	#	Purpose: initializes param with zeroes
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter
	#		manifold) together with imaginary component flagging 
	#		whether the hyper-parameter is fixed for purposes of estimation.
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Inputs:
	#		mdl: the specified sigex model, a list object
	#		data.ts: a T x N matrix ts object
	#	Outputs:
	#		list object with par.default and flag
	#		par.default: this is param, filled from a psi of all zeroes;
	#			see background.  Will have form specified by mdl.
	#		flag: string of ones, of length same as psi,
	#			with a 1 denoting that the corresponding hyper-parameter is to
	#			be estimated (by default, no parameters are fixed).
	#	Requires: sigex.zetalen, sigex.psi2par
	#
	####################################################################
 
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]
	A.mat <- matrix(0,N,N)
	A.mat[lower.tri(A.mat)] <- 1
	psi.len <- 0
	boundlist <- mdl[[5]]

	for(i in 1:length(mdl[[3]]))
	{	
		vrank <- mdl[[1]][[i]]
		D.dim <- length(vrank)
		L.dim <- sum(A.mat[,as.vector(vrank)])
		psi.len <- psi.len + D.dim + L.dim
		mdlType <- mdl[[2]][i]
		psi.len <- psi.len + sigex.zetalen(mdlType)
	}
	for(k in 1:N)
	{
		psi.len <- psi.len + dim(mdl[[4]][[k]])[2]
	}
	psi <- rep(0,psi.len) + 1i*rep(1,psi.len)
	par.default <- sigex.psi2par(psi,mdl,data.ts)
	flag <- Im(psi)

	return(list(par.default,flag))
}


