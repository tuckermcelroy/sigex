sigex.tstats <- function(mdl,psi,hess)
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
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter
	#		manifold) together with imaginary component flagging 
	#		whether the hyper-parameter is fixed for purposes of estimation.
	#     	The standard error has no division by T, because mvar.lik
	#		is T times the scale of the Whittle likelihood; multiply by 2 because 
	#		mvar.lik is -2 * log lik
	#	Inputs:
	#		mdl: the specified sigex model, a list object
	#		psi: see background 
	#		hess: Hessian matrix, which can be obtained from output of sigex.mlefit  
	#	Outputs:
 	#		tstats: vector of psi mles divided by standard error.  If a parameter
	#			is fixed during estimation, plus or minus infinity is returned.
	#			(This corresponds to a standard error of zero.)
	#
	####################################################################
  
	se <- rep(0,sum(Im(psi)))
	if(min(eigen(hess)$value) > 0) se <- sqrt(2*diag(solve(hess)))
	
	tstats <- Re(psi)
	tstats[Im(psi)==1] <- tstats[Im(psi)==1]/se
	tstats[Im(psi)==0] <- sign(tstats[Im(psi)==0])*Inf
	return(tstats)
}


