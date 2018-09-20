polymult <- function(a,b) 
{

	##########################################################################
	#
	#	polymult
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
	#	Purpose: compute the product of two polynomials
	#	Inputs:	
	#		a: vector of polynomial coefficients, where a[1] is the zeroth coefficient,
	#			a[2] is the first coefficient, etc.
	#		b: vector of polynomial coefficients, where b[1] is the zeroth coefficient,
	#			b[2] is the first coefficient, etc.
	#	Outputs: 
	#		c: vector of polynomial coefficients for c(z) = a(z)*b(z),
	#			where c[1] is the zeroth coefficient, c[2] is the first coefficient, etc.
	#
	##########################################################################################

	bb <- c(b,rep(0,length(a)-1))
	B <- toeplitz(bb)
	B[lower.tri(B)] <- 0
	aa <- rev(c(a,rep(0,length(b)-1)))
	prod <- B %*% matrix(aa,length(aa),1)
	return(rev(prod[,1]))
}