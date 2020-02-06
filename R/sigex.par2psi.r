sigex.par2psi <- function(param,mdl)
{

	##########################################################################
	#
	#	sigex.par2psi
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
	#	Purpose: transform param to psi
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Notes: this is a functional inverse to sigex.psi2par
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Inputs:
	#		param: see background.  Must have form specified by mdl
	#		mdl: the specified sigex model, a list object
	#	Outputs:
	#		psi: see background.
	#	Requires: sigex.par2zeta
	#
	####################################################################

	xi <- NULL
	zeta <- NULL
	N <- dim(as.matrix(param[[1]][[1]]))[1]
	for(i in 1:length(mdl[[3]]))
	{
		mdlType <- mdl[[2]][[i]]
		delta <- mdl[[3]][[i]]
		vrank <- mdl[[1]][[i]]
		L.mat <- param[[1]][[i]]
		D.mat <- param[[2]][[i]]
		new.xi <- c(L.mat[lower.tri(diag(N))[,as.vector(vrank)]],D.mat) 
		xi <- c(xi,new.xi)
		new.zeta <- sigex.par2zeta(param[[3]][[i]],mdlType)
		zeta <- c(zeta,new.zeta)
	}

	beta <- param[[4]]
	psi <- c(xi,zeta,beta)

	return(psi)
}
