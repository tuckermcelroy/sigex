sigex.tstats <- function(mdl,psi,hess,constraint)
{

	##########################################################################
	#
	#	sigex.tstats
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
	#	Purpose: computes t statistics for parameter estimates
	#	Background:
	#		param is the name for the model parameters entered into
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#   The standard error has no division by T, because mvar.lik
	#		is T times the scale of the Whittle likelihood; multiply by 2 because
	#		mvar.lik is -2 * log lik
	#	Inputs:
	#		mdl: the specified sigex model, a list object
	#		psi: see background
	#		hess: Hessian matrix, which can be obtained from output of sigex.mlefit
  #		constraint: matrix of the form [Q , C], with C (constraint.mat)
  #     the matrix of constraints and Q (constraint.vec) the vector
  #     of constraint constants, such that C psi = Q.
  #     Use NULL if there are no constraints
	#	Outputs:
 	#		tstats: vector of psi mles divided by standard error.  If a parameter
	#			constraints are used, this is taken into account.
	# Requires: sigex.eta2psi
  #
	####################################################################

	trans.mat <- diag(length(psi))
	if(length(constraint) > 0)
	{
	  constraint.mat <- constraint[,-1,drop=FALSE]
	  constraint.vec <- constraint[,1,drop=FALSE]
	  fixed.dim <- dim(constraint.mat)[1]
	  free.dim <- dim(constraint.mat)[2] - fixed.dim
	  trans.mat <- NULL
	  for(j in 1:free.dim)
	  {
	    eta <- diag(free.dim)[,j]
	    trans.mat <- cbind(trans.mat,sigex.eta2psi(eta,constraint))
	  }
	}

	se <- rep(0,length(psi))
	if(min(eigen(hess)$value) > 0)
	{
	  se <- sqrt(2*diag(trans.mat %*% solve(hess) %*% t(trans.mat)))
	}
	tstats <- psi/se
	return(tstats)
}


