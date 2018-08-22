sigex.blocktoep <- function(x.array)
{
	
	##########################################################################
	#
	#	sigex.blocktoep
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
	#	Purpose: generates lower triangular array for block Toeplitz matrix
	#	Background:	
	#		Given a sequence of matrix gamma(0), gamma(1),...,gamma(k)
	#		a block Toeplitz matrix has the form
	#		[ gamma(0), t(gamma(1)),...,t(gamma(k)) ]
	#		[ gamma(1), gamma(0),...,t(gamma(k-1)) ]
	#		[ gamma(k),...,gamma(0)]
	#	Notes: presumes that gamma(0) is symmetric
	#	Inputs:
	#		x.array: array of dimension N x N x H, where x.array[,,1] is gamma(0)
	#	Outputs:
	#		x.toep: array of dimension N x H x N x H, where x.toep[,j,,k]
	#			corresponds to gamma(j-k) if j >= k and is zero otherwise
	#
	####################################################################
  
	N <- dim(x.array)[1]
	H <- dim(x.array)[3]

	x.toep <- array(0,c(N,H,N,H))
	for(j in 1:H)
	{
		for(k in 1:H)
		{
			if(j >= k) { x.toep[,j,,k] <- x.array[,,j-k+1] }
		}
	}

	return(x.toep)
}
