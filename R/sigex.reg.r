sigex.reg <- function(mdl,series,reg)
{

	##########################################################################
	#
	#	sigex.reg
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
	#	Purpose: adds regressors to an existing model
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
	#	Notes: always use sigex.meaninit before calling this function.
	#		This function presumes that correct start and end dates for series
	#		of length T are known, and regressors are constructed accordingly!
	#		When regressors are entered, there is a check: the full differencing
	#		operator is applied, and if the magnitude of differenced regressors is small,
	#		they are deemed to be in the null space of the differencing operator,
	#		being non-identifiable -- hence they are omitted.
	#	Inputs:
	#		mdl: the specified sigex model, a list object.
	#			is your first component, then set mdl <- NULL
	#			mdl[[1]] is mdlK, gives ranks of white noise covariance matrix
	#			mdl[[2]] is mdlType, a list giving t.s. model class, order, and bounds
	#			mdl[[3]] is mdlDiff, gives delta differencing polynomials
	#		      mdl[[4]] is list of regressors by individual series
 	#		series: integer between 1 and N, the index of the individual series for
	#			which regressors are being added.
	#		reg: a one-column matrix of time series regressors, of length T.
  #     should have names attribute, as well as start and frequency
	#	Outputs:
	#		mdl: the updated sigex model, a list object
	#	Requires: sigex.delta
	#
	####################################################################

	mdlK <- mdl[[1]]
	mdlType <- mdl[[2]]
	mdlDiff <- mdl[[3]]
	mdlReg <- mdl[[4]]

	T <- dim(reg)[1]
	delta <-  sigex.delta(mdl,0)
	reg.diff <- filter(reg,delta,method="convolution",sides=1)[length(delta):T]
	if(sum(reg.diff^2) > 10^(-8))
	{
		mdlReg[[series]] <- ts(cbind(mdlReg[[series]][1:T,],reg[1:T,]),start=start(reg),
			frequency=frequency(reg),names=c(colnames(mdlReg[[series]]),colnames(reg)))
	}
	mdl <- list(ranks = mdlK,type = mdlType,diffop = mdlDiff,regress = mdlReg)

	return(mdl)
}

