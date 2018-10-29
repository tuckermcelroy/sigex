sigex.zetalen <- function(mdlType)
{

	##########################################################################
	#
	#	sigex.zetalen
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
	#	Purpose: computes the length of zeta
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
	#		mdlType: this is a component of mdl (the specified sigex model),
	#			 cited as mdl[[2]]
	#	Outputs:
	#		zetalen: the length of zeta
	#
	####################################################################

	mdlClass <- mdlType[[1]]
	mdlOrder <- mdlType[[2]]

	# ARMA  
	if(mdlClass %in% c("arma","arma.stab")) zeta.len <- sum(mdlOrder)

	# cycles
	if(mdlClass %in% c("bw","bw.stab","bal","bal.stab")) zetalen <- 2 

	return(zetalen)
}
