specFact <- function(poly)
{

	##########################################################################
	#
	#	specFact
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
	#	Purpose: compute the spectral factorization of a given p.d. sequence
	#	Background: any symmetric sequence has real-valued Fourier transform.
	#		If this is positive as well, the sequence is positive definite (p.d.),
	#		and can be factored as the magnitude squared of the Fourier transform
	#		of a causal power series, called the spectral factorization.  
	#		In the case that the p.d. sequence has finitely many terms, 
	#		it corresponds to the acf of an MA process and the MA polynomial
	#		is proportional to the spectral factorization.
	#	Inputs:	
	#		poly: a symmetric vector of coefficients, which is p.d.
	#	Outputs: 
	#		theta: the spectral factorization, such that
	#		poly(z,z^{-1}) = theta(z) * theta(z^{-1})
	#	Requires: polymult
	#
	##########################################################################################

	p <- length(poly)-1
	roots <- polyroot(poly)
	theta <- 1
	prod <- poly[p+1]
	toggle <- 1
	for(i in 1:p)
	{
		if (Mod(roots[i]) < 1)
		{
			theta <- polymult(theta,c(1,-roots[i]))
			prod <- prod/(-roots[i])
		} else {
		if (Mod(roots[i]) <= 1) 
		{
			if(Arg(roots[i]) < 0)
			{
				theta <- polymult(theta,c(1,-roots[i]))
				prod <- prod/(-roots[i])
			} else {
			if((Arg(roots[i]) == 0) && (toggle == 1))
			{	
				theta <- polymult(theta,c(1,-roots[i]))
				prod <- prod/(-roots[i])
				toggle <- -1*toggle
			} }
		} }
	}
	theta <- Re(theta)*sqrt(Re(prod))
	return(theta)
}
		  