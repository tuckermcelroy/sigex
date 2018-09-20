sigex.delta <- function(mdl,omits)
{

	##########################################################################
	#
	#	sigex.delta
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
	#	Purpose: compute differencing polynomial with factors omitted
	#	Background:	
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		The differencing polynomial for x_t is the product of
	#		the components' polynomials, so long as they are relatively prime
	#		(this is assumed).  Applying this product polynomial to x_t,
	#		the effect on a summand y_t is that it is differenced to 
	#		stationary w_t, but a remainder polynomial acts on w_t as well,
	#		given by the product of all other polynomials.
	#	Inputs:
	#		mdl: the specified sigex model, a list object.   
 	#		omits: indices of components that are to be omitted,
	#			when computing the product of differencing polynomials
	#			for the various components.
	#	Outputs:
	#		updated differencing polynomial,
	#		written in format c(delta0,delta1,...,deltad)
	#	Requires: polymult
	#
	####################################################################

	prod <- 1
	for(i in 1:length(mdl[[3]]))
	{
		polyn <- mdl[[3]][[i]]
		if (i %in% omits) polyn <- 1
		prod <- polymult(prod,polyn)		
	}
	return(prod)
}
