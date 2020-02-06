sigex.fixed <- function(data.ts,mdl,series,param,type)
{

	##########################################################################
	#
	#	sigex.fixed
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
	#	Purpose: computes all nontrend fixed regression effects 
	#	Background:
	#		x is a multivariate time series (N x T), and each individual series
	#		can have its distinct set of regressors.  So for each 1 <= j <= N,
	#		x[j,] is length T and has r_j number of length T regressors.
	#		There is a default regressor of polynomial time: suppose the
	#		time series has d unit roots (d >= 0), and this applies to each
	#		individual series (differencing polynomials are the same for all
	#		individual series in sigex).  Then the regressor t^d for 1 <= t <= T
	#		is the default "mean effect".  (Coefficients of lower order time
 	#		polynomial effects cannot be identified.)  When d=0, this is
	#		just the mean of the process.  (Although it need not be stationary
	#		when d=0, any other non-stationary latent components are assumed to
	#		have mean zero for identifiability.)  One can always add higher order
	#		time polynomial regressors, if desired.
	#		param is the name for the model parameters entered into 
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold) 
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		mdl: the specified sigex model, a list object. 
 	#		series: integer between 1 and N, the index of the individual series for 
	#			which regression effects are being computed.  
	#		param: see background.  Must have form specified by mdl
	#		type: a string designating the name of the regression effect
	#	Outputs:
	#		mean.mat: length T time series consisting of regression effects.
	#			This has the format X %*% beta, where each column of
	#			X is a regressor of "type", and beta consists of the corresponding
	#			regression parameters.
	#
	####################################################################
	
	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	ind <- 0
	mean.mat <- NULL
	for(k in 1:N)
	{
		reg.mat <- mdl[[4]][[k]]
		len <- dim(reg.mat)[2]
		if (k==series)
		{
			col.ind <- which(colnames(mdl[[4]][[k]])==type)
			if(length(col.ind) > 0) 
			{ 
				mean.mat <- as.matrix(reg.mat[,col.ind]) %*% 
					as.matrix(param[[4]][ind+col.ind]) 
				mean.mat <- ts(mean.mat,names=type)
			}
		}
		ind <- ind+len
	}
	
	return(mean.mat)
}
	



