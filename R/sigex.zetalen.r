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

	if(mdlType == "wn") { zetalen <- 0 }
	if(mdlType == "canonWN") { zetalen <- 0 }
	if(mdlType == "AR1") { zetalen <- 1 }
	if(mdlType == "MA1") { zetalen <- 1 }
	if(mdlType == "canonMA1") { zetalen <- 1 }
	if(mdlType %in% c("cycleBW1","cycleBW2","cycleBW3","cycleBW4","cycleBW5",
		"cycleBW6","cycleBW7","cycleBW8","cycleBW9","cycleBW10"))
		{ zetalen <- 2 }
	if(mdlType %in% c("canonCycleBW1","canonCycleBW2","canonCycleBW3",
		"canonCycleBW4","canonCycleBW5","canonCycleBW6","canonCycleBW7",
		"canonCycleBW8","canonCycleBW9","canonCycleBW10"))
		{ zetalen <- 2 }
	if(mdlType %in% c("cycleBAL1","cycleBAL2","cycleBAL3","cycleBAL4","cycleBAL5",
		"cycleBAL6","cycleBAL7","cycleBAL8","cycleBAL9","cycleBAL10"))
		{ zetalen <- 2 }
	if(mdlType %in% c("canonCycleBAL1","canonCycleBAL2","canonCycleBAL3",
		"canonCycleBAL4","canonCycleBAL5","canonCycleBAL6","canonCycleBAL7",
		"canonCycleBAL8","canonCycleBAL9","canonCycleBAL10"))
		{ zetalen <- 2 }
#	if(mdlType == "VAR1") { zetalen <- N^2 }
	if(mdlType == "ARMA22") { zetalen <- 4 }
	if(mdlType == "damped") { zetalen <- 1 + as.integer(N*(N+1)/2) }

	return(zetalen)
}
