sigex.canonize <- function(ma.coef,ar.coef)
{

	##########################################################################
	#
	#	sigex.canonize
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
	#	Purpose: compute canonization of a given ARMA process
	#	Background:	
	#		An ARMA or ARIMA process X_t has form
	#			X_t = theta(B)/phi(B) eps_t
	#		where eps_t is white noise, and we allow phi(B) to have unit roots.
	#		Canonization seeks a new theta*(B) such that the spectral density
	#		corresponding to theta*(B)/phi(B) is non-invertible (i.e. has a zero).
	#		The minimum value of the original spectrum is subtracted off to
	#		get the canonized spectrum.	
	#	Inputs:
	#		ma.coef: q coefficients of MA polynomial with unit constant coefficient
	#		ar.coef: p coefficients (minus convention) of AR polynomial with unit constant coefficient
	#	Outputs:
	#		ma.stab: coefficients of MA polynomial (without a unit constant coefficient),
	#			 the degree is max(p,q)
	#	Requires: polymult, polysum, specFact
	#
	####################################################################
	
	q <- length(ma.coef)
	p <- length(ar.coef)
	r <- 2*(q+p)
	thresh <- 10^(-8)

	# find candidate frequencies of critical points
	phi.poly <- c(1,-1*ar.coef)
	phi.grad.poly <- phi.poly*seq(0,p)
	phi.rev.poly <- rev(phi.poly)
	phi.grad.rev.poly <- rev(phi.grad.poly)
	theta.poly <- c(1,ma.coef)
	theta.grad.poly <- theta.poly*seq(0,q)
	theta.rev.poly <- rev(theta.poly)
	theta.grad.rev.poly <- rev(theta.grad.poly)
	numer.phi <- polymult(phi.grad.poly,phi.rev.poly) - polymult(phi.poly,phi.grad.rev.poly)
	numer.phi <- polymult(numer.phi,polymult(theta.poly,theta.rev.poly))
	numer.theta <- polymult(theta.grad.poly,theta.rev.poly) - polymult(theta.poly,theta.grad.rev.poly)
	numer.theta <- polymult(numer.theta,polymult(phi.poly,phi.rev.poly))
	numer <- -1*numer.phi + numer.theta
	my.roots <- polyroot(numer)
	lambdas <- c(0,pi)
	for(i in 1:r)
	{
		if( abs(Mod(my.roots[i])-1) < thresh ) { lambdas <- c(lambdas, Arg(my.roots[i])) }
	}

	# determine global minimum
	min.val <- Inf
	min.ind <- 0
	for(j in 1:length(lambdas))
	{
		new.val <- Mod(sum(exp(-1i*lambdas[j]*seq(0,q))*theta.poly))^2/Mod(sum(exp(-1i*lambdas[j]*seq(0,p))*phi.poly))^2
		if(new.val < min.val) { min.val <- new.val; min.ind <- j }
	}

	# compute new MA polynomial
	numer.acf <- polysum(polymult(theta.poly,rev(theta.poly)),-1*min.val*polymult(phi.poly,rev(phi.poly)))
 	ma.stab <- Re(specFact(numer.acf))	

	return(ma.stab)
}


