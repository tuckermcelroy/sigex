sigex.seascheck <- function(signal,period,tolerance,diffs)
{

	##########################################################################
	#
	#	sigex.seascheck
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
	#	Purpose: computes a seasonality diagnostic
	#	Background:	
	#		For each series in the estimated signal,
	#		there is an OLS AR fit (the signal gets differenced
 	#		(1-B)^diffs times).  The the roots are computed, and
	#		the phase is compared to seasonal roots, with magnitudes
	#		close to unity indicating "visual significance".  
 	#		The tolerance level determines how close the phase should
	#		be for the seasonality.  Heuristically, if the root magnitude
	#		is less than 1.03, then there is "visual significance", and 
	#		hence seasonality may be present.  
	#		If no root of a particular seasonal frequency
	#		is found, the function returns Inf.  
	#		Returns floor(period/2) values for each series 
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
	#		Signal extraction diagnostics are computed, which are based
	#		on autocorrelations of differenced estimated signals, compared
	#		to their expected variability.  Blakely and McElroy (2017 Econometric Reviews)
	#		discusses the theory in generality; this version does
	#		not take parameter uncertainty into account.
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter
	#		manifold) together with imaginary component flagging 
	#		whether the hyper-parameter is fixed for purposes of estimation.
	#	Inputs:
	#		signal: a T x N matrix ts object of the extracted signal
	#		period: the number of seasons per year
	#		tolerance: see background
	#		diffs: integer number of regular differences to be applied.
	#	Outputs:
	#		mag.mat: matrix with N columns and number of rows given
 	#			by period/2, giving the magnitude for the given series
	#			at each seasonal root, so long as the phase is within
	#			the tolerance level.
	#
	####################################################################
 
	x <- t(signal)
	N <- dim(x)[1]
	T <- dim(x)[2]

	mag.mat <- NULL
	for(i in 1:N) 
	{
		mags <- NULL
		diffed.x <- x[i,]
		if(diffs > 0) { diffed.x <- diff(x[i,],differences = diffs) }
		my.arfit <- ar.ols(diffed.x)
		my.roots <- polyroot(c(1,-1*my.arfit$ar[,,1]))
		polar.roots <- cbind(round(Mod(my.roots),digits=3),
			signif(2*pi/Arg(my.roots),digits=6),round(period*Arg(my.roots)/(2*pi),digits=2))
		for(k in 1:(period/2))
		{
			mag <- polar.roots[abs(polar.roots[,3] - k) < tolerance,1]
			if(length(mag)==0) mag <- Inf
			if(length(mag)>1) mag <- min(mag)
			mags <- c(mags,mag)
		}
		mag.mat <- cbind(mag.mat,mags)
	}
				
	return(mag.mat)

}
