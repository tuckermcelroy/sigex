sigex.lpwk <- function(data.ts,param,mdl,trendcyclecomp,grid,len,cutoff,trunc)
{

	##########################################################################
	#
	#	sigex.lpwk
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

	################# Documentation ############################################
	#
	#	Purpose: computes signal extraction filter coefficients
	#		for trend and cycle, by combining WK filter for trend-cycle
	#		(specified by trendcyclecomp) with LP filter of cutoff.
	#	Background:	
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		We have a model for each w_t process, and can compute its
	#		autocovariance function (acf), and denote its autocovariance
	#		generating function (acgf) via gamma_w (B).
	#		The signal extraction filter for y_t is determined from
	#		this acgf and delta.  
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Notes: take grid >> len, else numerical issues arise
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background.  Must have form specified by mdl
	#		mdl: the specified sigex model, a list object
	#		trendcyclecomp: is the (single) index of the trend-cycle component
	#		sigcomps: provides indices of a desired component that
	#			is disjoint from trend-cycle, so that MSEs of 
	#			trend+sigcomps and cycle+sigcomps are computed.
	#		 	(Pass in sigcomps = NULL to just get trend and cycle MSEs.)
	#		grid: desired number of frequencies for spectrum calculations
	#		len: max index of the filter coefficients
	#		cutoff: is a number between 0 and pi, with all frequencies < cutoff preserved
	#		trunc: truncation index for LP filter
	#	Outputs:
	#		list object with psi.lpwk and psi.ilpwk
	#		psi.lpwk: array of dimension c(N,N,2*trunc+1), filter coefficients for trend
	#		psi.ilpwk: array of dimension c(N,N,2*trunc+1), filter coefficients for cycle
	#	Requires: sigex.wk
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	target <- array(diag(N),c(N,N,1))
	wkcoeff <- sigex.wk(data.ts,param,mdl,trendcyclecomp,target,FALSE,grid,len)
	mid <- trunc+len
	lpcoeff <- c(cutoff/pi,sin(cutoff*seq(1,mid))/(pi*seq(1,mid)))	
	lpcoeff <- c(rev(lpcoeff),lpcoeff[-1]) 
	ilpcoeff <- c(1-cutoff/pi,-sin(cutoff*seq(1,mid))/(pi*seq(1,mid)))	
	ilpcoeff <- c(rev(ilpcoeff),ilpcoeff[-1]) 

	psi.lpwk <- matrix(wkcoeff[[1]],nrow=N) %*% (lpcoeff[(mid+1+len):(mid+1-len)] %x% diag(N))
	psi.ilpwk <- matrix(wkcoeff[[1]],nrow=N) %*% (ilpcoeff[(mid+1+len):(mid+1-len)] %x% diag(N))
	for(k in 1:trunc)
	{
		next.coeff <- matrix(wkcoeff[[1]],nrow=N) %*% (lpcoeff[(mid+1+k+len):(mid+1+k-len)] %x% diag(N))
		psi.lpwk <- cbind(psi.lpwk,next.coeff)
		next.coeff <- matrix(wkcoeff[[1]],nrow=N) %*% (ilpcoeff[(mid+1+k+len):(mid+1+k-len)] %x% diag(N))
		psi.ilpwk <- cbind(psi.ilpwk,next.coeff)
		prev.coeff <- matrix(wkcoeff[[1]],nrow=N) %*% (lpcoeff[(mid+1-k+len):(mid+1-k-len)] %x% diag(N))
 		psi.lpwk <- cbind(prev.coeff,psi.lpwk)
		prev.coeff <- matrix(wkcoeff[[1]],nrow=N) %*% (ilpcoeff[(mid+1-k+len):(mid+1-k-len)] %x% diag(N))
 		psi.ilpwk <- cbind(prev.coeff,psi.ilpwk)
	}

	psi.lpwk <- array(psi.lpwk,c(N,N,(2*trunc+1)))
	psi.ilpwk <- array(psi.ilpwk,c(N,N,(2*trunc+1)))

	return(list(psi.lpwk,psi.ilpwk))
}


