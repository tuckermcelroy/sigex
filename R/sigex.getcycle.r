sigex.getcycle <- function(cycle.order,rho,omega)
{

	##########################################################################
	#
	#	sigex.getcycle
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
	#	Purpose: compute AR and MA polynomials for Butterworth and Balanced cycles
	#	Background:	
	#		A stochastic cycle is a particular parametrization of an ARMA process,
	#		governed by its order, its persistency (rho), and its period (2/omega)	
	#	Inputs:
	#		cycle.order: integer order of the cycle
	#		rho: persistency, a value between 0 and 1
	#		omega: frequency in units of pi, a value between 0 and 1.
	#			(The period is 2/omega)
	#	Outputs:
	#		list object with cycle.AR and cycle.MA
	#		cycle.AR: full AR polynomial in format c(1,-phi1,...,-phip)
	#		cycle.MA: full MA polynomial in format c(1,theta1,...,thetaq)
	#	Requires: polymult
	#
	####################################################################
	
	cycle.AR <- 1
	cycle.MA <- 1
	if(cycle.order > 0)
	{
		if(omega==0) {
		for(l in 1:cycle.order)
		{
			cycle.AR <- polymult(cycle.AR,c(1,-rho))
		} }
		if(omega==1) {
		for(l in 1:cycle.order)
		{
			cycle.AR <- polymult(cycle.AR,c(1,rho))
		} }
		if((omega > 0) & (omega < 1)) {
		for(l in 1:cycle.order)
		{
			cycle.AR <- polymult(cycle.AR,c(1,-2*rho*cos(pi*omega),rho^2))
			cycle.MA <- polymult(cycle.MA,c(1,-1*rho*cos(pi*omega)))
		} 
		#cycle.acf <- polymult(cycle.MA,rev(cycle.MA))
		}
	}

	return(list(cycle.AR,cycle.MA))
}


