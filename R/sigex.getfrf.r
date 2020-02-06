sigex.getfrf <- function(data.ts,param,mdl,sigcomps,plotit=TRUE,grid)
{

	##########################################################################
	#
	#	sigex.getfrf
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
	#	Purpose: computes frf for desired signal and plots it
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
	#		sigcomps: indices of the latent components composing the signal
	#		plotit: boolean flag for whether frf should be plotted;
	#			only plots if N <= 3
	#		grid: desired number of frequencies for spectrum calculations
	#	Outputs:
	#		frf.comp:  array of dimension c(N,N,grid), with complex number entries 
	#	Requires: sigex.frf
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]

	frf.comp <- sigex.frf(data.ts,param,mdl,sigcomps,grid)

	if(plotit && (N <= 4)) {
	par(mfrow = c(N,N))
	for(i in 1:N) 
	{
		for(j in 1:N)
		{
			plot(ts(Re(frf.comp[i,j,]),start=0,frequency=grid),
				xlab="Cycles",ylab="")
		}
	} }

	return(frf.comp)
}


