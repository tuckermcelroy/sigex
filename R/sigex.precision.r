sigex.precision <- function(data.ts,param,mdl,sigcomps)
{

	##########################################################################
	#
	#	sigex.precision
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
	#	Purpose: compares signal extraction MSE arising from
	#	 	multivariate fit with the implied univariate fit.
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
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background
	#		mdl: the specified sigex model, a list object
	#		sigcomps: indices of the latent components composing the signal
	#	Outputs:
	#		ratio of mse to mse.uni, the signal extraction MSE for 
	#			multivariate and univariate models
	#	Requires: sigex.mvar2uvar, sigex.signal
	#
	####################################################################
 
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	out <- sigex.mvar2uvar(data.ts,param,mdl)
	mdl.uni <- out[[1]]
	par.uni <- out[[2]]

	signal <- sigex.signal(data.ts,param,mdl,sigcomps)
	signal.uni <- sigex.signal(data.ts,par.uni,mdl.uni,sigcomps)

	mse <- matrix(diag(signal[[2]]),T,N)
	mse.uni <- matrix(diag(signal.uni[[2]]),T,N)

	return(mse/mse.uni)
}
