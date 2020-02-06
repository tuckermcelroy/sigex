sigex.mvar2uvar <- function(data.ts,param,mdl)
{

	##########################################################################
	#
	#	sigex.mvar2uvar
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
	#	Purpose: transform multivariate parameters to implied univariate form
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background
	#		mdl: the specified sigex model, a list object
	#	Outputs:
	#		mdl.uni: output univariate model
	#		univ.param: param for univariate model
	#	Requires: sigex.default
	#
	####################################################################
 
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	mdl.uni <- mdl
	for(i in 1:length(mdl.uni[[1]])) { mdl.uni[[1]][[i]] <- seq(1,N) }
	default.param <- sigex.default(mdl.uni,x)
	univ.param <- param
	for(i in 1:length(mdl[[3]]))
	{
		L.psi <- param[[1]][[i]]
		D.psi <- param[[2]][[i]]
		k.order <- length(D.psi)
	 	S.mat <- L.psi %*% diag(exp(D.psi),nrow=k.order) %*% t(L.psi)
		univ.param[[2]][[i]] <- log(diag(S.mat))				
		univ.param[[1]][[i]] <- default.param[[1]][[i]]
	}

	return(list(mdl.uni,univ.param))
}
