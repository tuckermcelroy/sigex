sigex.meaninit <- function(mdl,data.ts,d)
{

	##########################################################################
	#
	#	sigex.meaninit
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
	#	Purpose: adds trend regressors to an existing model
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
	#	Notes: always use this function when setting up the model.
	#		First make all calls to sigex.add, then call sigex.meaninit,
	#		and then add additional regressors (if needed) with sigex.reg.   
	#	Inputs:
	#		mdl: the specified sigex model, a list object. 
	#			mdl[[1]] is mdlK, gives ranks of white noise covariance matrix
	#			mdl[[2]] is mdlType, a list giving t.s. model class, order, and bounds
	#			mdl[[3]] is mdlDiff, gives delta differencing polynomials
	#		      mdl[[4]] is list of regressors by individual series
	#		data.ts: a T x N matrix ts object
	#		d: order of time polynomial trend regressor desired, labeled as "Trend".
	#			However, if a trend component exists in the model,
	#			then this d is ignored, and the number of unit roots
	#			is instead used to determine d.
	#	Outputs:
	#		mdl: the updated sigex model, a list object
	#	Requires: sigex.whichtrend, sigex.delta
	#
	####################################################################

	mdlK <- mdl[[1]]
	mdlType <- mdl[[2]]
	mdlDiff <- mdl[[3]]
	mdlReg <- mdl[[4]]

	if(length(sigex.whichtrend(mdl))==1) 
	{
		delta <- mdl[[3]][[sigex.whichtrend(mdl)]]
		d <- length(delta) - 1 - sum(abs(Arg(polyroot(delta))) > 10^(-8))
	}
	x <- t(data.ts)
	T <- dim(x)[2]
	N <- dim(x)[1]
	delta <- sigex.delta(mdl,0)
	for(series in 1:N)
	{ 
		mdlReg[[series]] <- ts(as.matrix(seq(1,T)^d),start=start(data.ts),
			frequency=frequency(data.ts),names="Trend") 
		d.inc <- d
		while(d.inc > 0)
		{
			d.inc <- d.inc-1
			reg <- ts(as.matrix(seq(1,T)^d.inc),start=start(data.ts),
				frequency=frequency(data.ts),names="Trend")
			reg.diff <- filter(reg,delta,method="convolution",
				sides=1)[length(delta):T] 
			if(sum(reg.diff^2) > 10^(-8)) 
			{
				mdlReg[[series]] <- ts(cbind(mdlReg[[series]],reg),
					start=start(reg),frequency=frequency(reg),
					names=c(colnames(mdlReg[[series]]),colnames(reg)))
			}
		}
	}
	mdl <- list(ranks = mdlK,type = mdlType,diffop = mdlDiff,regress = mdlReg)
	return(mdl)
}

