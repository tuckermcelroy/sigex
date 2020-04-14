sigex.glr <- function(data.ts,psi.nested,psi.nesting,mdl.nested,mdl.nesting)
{

	##########################################################################
	#
	#	sigex.glr
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
	#	Purpose: computes the difference of -2*log(Gaussian likelihood) for
	#		two models, the nested lik minus nesting lik
	#	Background:	
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Format: psi has three portions, psi = [xi,zeta,beta]
	#		xi ~ all hyper-parameters for covariance matrices
	#		zeta ~ all hyper-parameters for t.s. models
	#		beta ~ all regression parameters
	#	Notes: handles missing values in data.ts, which are indicated by 1i.
	#		This test can be applied to non-nested models,
	#		but the distribution won't be chi^2
	#	Inputs:
	#		data.ts: a T x N matrix ts object; any missing values 
	#			must be encoded with 1i in that entry
	#		psi.nested: see background; psi of the nested model
	#		psi.nesting: see background; psi of the nesting model
	#		mdl.nested: the specified nested sigex model, a list object
	#		mdl.nesting: the specified nesting sigex model, a list object
	#	Outputs:
	#		glr: difference of log likelihoods
	#		dof: degrees of freedom, given by different in number of parameters
	#	Requires: sigex.lik
	#
	####################################################################
 
  debug <- FALSE
	glr <- sigex.lik(psi.nested,mdl.nested,data.ts,debug) - 
		sigex.lik(psi.nesting,mdl.nesting,data.ts,debug)
	dof <- length(psi.nesting) - length(psi.nested)

	return(c(glr,dof))
}
