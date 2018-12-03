sigex.add <- function(mdl,vrank,class,order,bounds,name,delta)
{

	##########################################################################
	#
	#	sigex.add
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
	#	Purpose: build the model by adding on another latent component
	#	Background:	
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		We have a model for each w_t process, which is specified
	#		through the ranks (indices of non-zero Schur complements,
	#		cf. background for sigex.param2gcd) of the white noise 
	#		covariance matrix; also there is the model type, which 
	#		denotes the specification of the t.s. model for w_t;
	#		all the regressors, which are specified by individual time series
	#		rather than by latent component, and must have length T;
	#		pre-specified bounds for cyclical parameters, for each component,
	#		if applicable.
	#	Inputs:
	#		mdl: the specified sigex model, a list object.  This is whatever
	#			model you have so far specified, to which you will be adding
	#			model structure corresponding to the new component.  If this
	#			is your first component, then set mdl <- NULL
	#			mdl[[1]] is mdlK, gives ranks of white noise covariance matrix
	#			mdl[[2]] is mdlType, a list giving t.s. model class, order, and bounds
	#			mdl[[3]] is mdlDiff, gives delta differencing polynomials
	#		      mdl[[4]] is list of regressors by individual series
	#		vrank: vector of integers between 1 and N, corresponding
	#			to indices of non-zero Schur complements in the GCD
	#			of the innovations' covariance matrix for the new latent component
	#		class: character string of t.s. model type for the new latent component
	#		order: vector of model order
	#	      bounds: four numbers, gives bounds for rho and omega, 
	#			the cycle parameters of the new latent component
	#			rho lies in (bounds[1],bounds[2])
	#			omega lies in (bounds[3],bounds[4])
	#		name: character string giving the latent component's name
	#		delta: differencing polynomial (corresponds to delta(B) in Background)
	#			written in format c(delta0,delta1,...,deltad)
	#	Outputs:
	#		mdl: the updated sigex model, a list object
	#
	####################################################################

	mdlK <- mdl[[1]]
	mdlType <- mdl[[2]]
	mdlDiff <- mdl[[3]]
	mdlReg <- mdl[[4]]

	rank.null <- NULL
	if(length(vrank)==1) rank.null <- 0
	mdlK[[length(mdlK)+1]] <- c(rank.null,vrank)
	if(length(vrank)==1) mdlK[[length(mdlK)]] <- mdlK[[length(mdlK)]][-1]
	mdlType[[length(mdlType)+1]] <- list(class,order,bounds,name)
	delta.null <- NULL
	if(length(delta)==1) delta.null <- 0
	mdlDiff[[length(mdlDiff)+1]] <- c(delta.null,delta)
	if(length(delta)==1) mdlDiff[[length(mdlDiff)]] <- mdlDiff[[length(mdlDiff)]][-1]

	mdl <- list(ranks = mdlK,type = mdlType,diffop = mdlDiff,regress = mdlReg)
	return(mdl)
}

