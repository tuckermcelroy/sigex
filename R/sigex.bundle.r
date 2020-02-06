sigex.bundle <- function(data.ts,transform,mdl,psi)
{

	##########################################################################
	#
	#	sigex.bundle
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
	#	Purpose: bundles into a list object several features
	#	Background:
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		data.ts: a T x N matrix ts object, 
	#			corresponding to N time series of length T
	#		transform: a character indicating an instantaneous 
	#			transformation to be applied; current options are
	#			"none", "log", and "logistic"
 	#		mdl: the specified sigex model, a list object
	#		psi: see background.  
	#	Outputs:
	#		analysis: a list object of the inputs
	#
	####################################################################

	analysis <- list(data = data.ts,transform = transform,model = mdl,psi = psi)

	return(analysis)
}

