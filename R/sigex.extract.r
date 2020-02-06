sigex.extract <- function(data.ts,filter,mdl,param)
{

	##########################################################################
	#
	#	sigex.extract
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
	#	Purpose: computes signal extraction estimates with two standard errors
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
	#		filter: list object corresponding to the output of sigex.signal,
	#			a list object of f.mat and v.mat.
	#			f.mat: array of dimension c(T,N,T,N), where f.mat[,j,,k]
	#				is the signal extraction matrix that utilizes input series k
	#				to generate the signal estimate for series j.
	#			v.mat: array of dimension c(T,N,T,N), where v.mat[,j,,k]
	#				is the error covariance matrix arising from input series k
	#				used to generate the signal estimate for series j.
	#		mdl: the specified sigex model, a list object
	#		param: see background.  Must have form specified by mdl
	#	Outputs:
	#		list object with extract, upp, and low
	#		extract: T x N matrix of the signal estimates
	#		upp: as extract, plus twice the standard error
	#		low: as extract, minus twice the standard error
	#
	####################################################################
 
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	# subtract regression effects
	ind <- 0
	data.diff <- data.ts
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		data.diff[,k] <- data.diff[,k] - reg.mat %*% 
			as.matrix(param[[4]][(ind+1):(ind+len)])
		ind <- ind+len
	}
	xvec <- matrix(t(data.diff),nrow=N*T,ncol=1)
	 
  extract <- filter[[1]] %*% xvec
	extract <- t(matrix(extract,nrow=N,ncol=T))

	mse <- t(matrix(diag(filter[[2]]),N,T))
	upp <- extract + 2*sqrt(mse)
	low <- extract - 2*sqrt(mse)
	 
	return(list(extract,upp,low))
}
