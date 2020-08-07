sigex.reduce <- function(data.ts,param,mdl,thresh,modelflag)
{

	##########################################################################
	#
	#	sigex.reduce
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
	#	Purpose: determine a reduced rank model from a given fitted model
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
	#		param is the name for the model parameters entered into
	#		a list object with a more intuitive structure, whereas
	#		psi refers to a vector of real numbers containing all
	#		hyper-parameters (i.e., reals mapped bijectively to the parameter	manifold)
	#	Inputs:
	#		data.ts: a T x N matrix ts object
	#		param: see background
	#		mdl: the specified sigex model, a list object.  This is whatever
	#			model you have so far specified, to which you will be adding
	#			model structure corresponding to the new component.  If this
	#			is your first component, then set mdl <- NULL
	#			mdl[[1]] is mdlK, gives ranks of white noise covariance matrix
	#			mdl[[2]] is mdlType, a list giving t.s. model class, order, and bounds
	#			mdl[[3]] is mdlDiff, gives delta differencing polynomials
	#		      mdl[[4]] is list of regressors by individual series
	#		thresh: lower bound on Schur complements
	#		modelFlag: when TRUE, small Schur complements imply rank reduction
	#			in the new model.  When modelFlag is FALSE, small Schur
	#			complements are replaced by exp(thresh)
	#	Outputs:
	#		mdl.red: the new sigex model, a list object
	#		par.red: the new param for the new model
	#	Requires: sigex.par2psi, sigex.conditions, sigex.add,
	#		sigex.meaninit, sigex.renderpd
	#
	####################################################################

	x <- t(data.ts)
	N <- dim(x)[1]
	T <- dim(x)[2]

	psi.red <- sigex.par2psi(param,mdl)
	log.conds <- log(sigex.conditions(data.ts,psi.red,mdl))

	if(modelflag) {
	par.red <- param
	mdl.red <- NULL
	for(j in 1:length(mdl[[3]])) {
		ranks <- seq(1,N)[log.conds[j,] > thresh]
		mdlType <- mdl[[2]][[j]]
		mdl.red <- sigex.add(mdl.red,ranks,mdlType[[1]],mdlType[[2]],
			mdlType[[3]],mdlType[[4]],mdl[[3]][[j]])
		par.red[[1]][[j]] <- as.matrix(param[[1]][[j]][,ranks])
		par.red[[2]][[j]] <- param[[2]][[j]][ranks]
	}
	mdl.red <- sigex.meaninit(mdl.red,data.ts,0)
	mdl.red[[4]] <- mdl[[4]] } else {
	par.red <- param
	mdl.red <- mdl
	for(j in 1:length(mdl[[3]]))
	{
		L.mat <- par.red[[1]][[j]]
		D.psi <- par.red[[2]][[j]]
		D.new <- sigex.renderpd(L.mat,D.psi,thresh)
		par.red[[2]][[j]] <- D.new
	} }

	return(list(mdl.red,par.red))
}

