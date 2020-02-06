sigex.lpfiltering <- function(mdl,data.ts,trendcyclecomp,sigcomps,psi,cutoff,grid,window,trunc,trendFlag)
{

	##########################################################################
	#
	#	sigex.lpfiltering
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
	#	Purpose: computes signal extraction estimates with uncertainty
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
	#		this acgf and delta.  The error spectral density calculations
	#		are found in: 
	#		"Casting Vector Time Series: Algorithms for Forecasting,
	#		Imputation, and Signal Extraction," McElroy (2018).
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Notes:
	#		Starts with LP an ideal low-pass filter,
	#		and applies LP*WK filter with cutoff parameter  to
	#		each component, where the filter has been truncated to 
	#		length 2*window + 1.  The output will be LP*WK filter applied to data,
	#	  	which is forecast and aftcast extended by window units, covering time points
 	#	  	1-window ..., T+window.  This gives trend and cycle at times 1,...,T.
	#	      Take grid >> window, else numerical issues arise
	#	Inputs:
	#		mdl: the specified sigex model, a list object
	#		data.ts: a T x N matrix ts object
	#		trendcyclecomp: is the (single) index of the trend-cycle component
	#		sigcomps: provides indices of a desired component that
	#			is disjoint from trend-cycle, so that MSEs of 
	#			trend+sigcomps and cycle+sigcomps are computed.
	#		 	(Pass in sigcomps = NULL to just get trend and cycle MSEs.)
	#		psi: see background. 
	#		cutoff: is a number between 0 and pi, with all frequencies < cutoff preserved
	#		grid: desired number of frequencies for spectrum calculations
	#		window: max index of the filter coefficients
	#		trunc: truncation index for LP filter
	#		trendFlag: boolean flag, TRUE for trend+signal, else get cycle+signal
	#	Outputs:
	#		list object with lp.signal, upp, and low
	#		lp.signal: T x N matrix of the signal estimates
	#		upp: as lp.signal, plus twice the standard error
	#		low: as lp.signal, minus twice the standard error
	#	Requires: sigex.psi2par, sigex.lpwk, sigex.lpmse, sigex.cast, sigex.wkextract2
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]
	param <- sigex.psi2par(psi,mdl,data.ts)

	lpwk.out <- sigex.lpwk(data.ts,param,mdl,trendcyclecomp,grid,window,cutoff,trunc)
	lpwk.filter <- lpwk.out[[1]]
	ilpwk.filter <- lpwk.out[[2]]

	lp.mse <- sigex.lpmse(param,mdl,trendcyclecomp,sigcomps,grid,cutoff)
	if(trendFlag) 
	{ 
		lp.mse <- lp.mse[[1]] 
		psi.filter <- lpwk.filter
	} else 
	{ 
		lp.mse <- lp.mse[[2]] 
		psi.filter <- ilpwk.filter
	}	

	if(length(sigcomps) > 0)
	{
		extract.signal <- sigex.wkextract2(psi,mdl,data.ts,sigcomps,grid,window,0,NULL,FALSE)
	}
	if((trunc) > 0) { 
		leads <- c(-rev(seq(0,trunc-1)),seq(1,T),seq(T+1,T+trunc))
	} else { leads <- seq(1,T) }
	data.ext <- t(sigex.cast(psi,mdl,data.ts,leads,FALSE))

	lp.signal <- NULL
	upp <- NULL
	low <- NULL
	for(j in 1:N) 
	{
		output.j <- rep(0,T)
		for(k in 1:N)
		{
			output.k <- filter(data.ext[,k],psi.filter[j,k,],
				method="convolution",sides=2)
			output.k <- output.k[(trunc+1):(trunc+T)]
			output.j <- output.j + output.k
		}
		if(length(sigcomps) > 0) { output.j <- output.j + extract.signal[[1]][,j] }
		lp.signal <- cbind(lp.signal,output.j)
	 	mse <- lp.mse[j,j]
	 	upp <- cbind(upp,output.j + 2*sqrt(mse))
		low <- cbind(low,output.j - 2*sqrt(mse))
	}	
	
	return(list(lp.signal,upp,low))
}



