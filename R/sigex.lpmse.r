sigex.lpmse <- function(param,mdl,trendcyclecomp,sigcomps,grid,cutoff)
{

	##########################################################################
	#
	#	sigex.lpmse
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
	#	Purpose: computes signal extraction mse arising from LP filtering
	#		of trend-cycle, with low pass cutoff, for trend and cycle each.
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
	#		this acgf and delta.  The error spectral density calculations
	#		are found in: 
	#		"Casting Vector Time Series: Algorithms for Forecasting,
	#		Imputation, and Signal Extraction," McElroy (2018).
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		param: see background.  Must have form specified by mdl
	#		mdl: the specified sigex model, a list object
	#		trendcyclecomp: is the (single) index of the trend-cycle component
	#		sigcomps: provides indices of a desired component that
	#			is disjoint from trend-cycle, so that MSEs of 
	#			trend+sigcomps and cycle+sigcomps are computed.
	#		 	(Pass in sigcomps = NULL to just get trend and cycle MSEs.)
	#		grid: desired number of frequencies for spectrum calculations
	#		cutoff: is a number between 0 and pi, with all frequencies < cutoff preserved
	#	Outputs:
	#		list object with mse.trend and mse.cycle
	#		mse.trend: N x N matrix, MSE of trend
	#		mse.cycle: N x N matrix, MSE of cycle
	#	Requires: sigex.spectra, sigex.delta, sigex.wkmse
	#
	####################################################################

quad <- function(z)
{
	# quadrature of z 
	len <- length(z)
	out <- (z[1]+z[len])/2 + sum(z[2:(len-1)])
	return(out/len)
}

	N <- dim(as.matrix(param[[1]][[1]]))[1]
	lambda <- pi*seq(0,grid)/grid
	f.sig <- t(rep(0,grid+1) %x% diag(N))
	for(i in 1:length(mdl[[3]]))
	{
		L.par <- param[[1]][[i]]
		D.par <- param[[2]][[i]]
		delta <- sigex.delta(mdl,i)
		f.comp <- sigex.spectra(L.par,D.par,mdl,i,param[[3]][[i]],delta,grid)
		if(i %in% trendcyclecomp) f.sig <- f.sig + matrix(f.comp,nrow=N)
	}
	frf.wksig <- sigex.wkmse(data,param,mdl,trendcyclecomp,grid)
	frf.wksig <- matrix(frf.wksig,nrow=N)

	# get terms of MSE 
	allcomps <- seq(1,length(mdl[[3]]))
 	seminoisecomps <- allcomps[!allcomps %in% c(trendcyclecomp,sigcomps)]
	if(length(sigcomps)==0) { frf.wksig <- matrix(0,nrow=N,ncol=(N*(grid+1))) } else {
		frf.wksig <- sigex.wkmse(data,param,mdl,sigcomps,grid)
		frf.wksig <- matrix(frf.wksig,nrow=N) }
	if(length(seminoisecomps)==0) { frf.wknoise <- matrix(0,nrow=N,ncol=N*(grid+1)) } else {
		frf.wknoise <- sigex.wkmse(data,param,mdl,seminoisecomps,grid)
		frf.wknoise <- matrix(frf.wknoise,nrow=N) }
	integrand1 <- rep(0,length(lambda))
	integrand1[lambda < cutoff] <- 1
	integrand2 <- 1 - integrand1

	integrand <- t(integrand2 %x% matrix(1,N,N)) * frf.wksig +
				t(integrand1 %x% matrix(1,N,N)) * frf.wknoise
	integrand <- array(integrand,c(N,N,(grid+1)))
	mse.trend <- Re(apply(integrand,c(1,2),quad))

	integrand <- t(integrand1 %x% matrix(1,N,N)) * frf.wksig +
				t(integrand2 %x% matrix(1,N,N)) * frf.wknoise
	integrand <- array(integrand,c(N,N,(grid+1)))
	mse.cycle <- Re(apply(integrand,c(1,2),quad))

	return(list(mse.trend,mse.cycle))
}
