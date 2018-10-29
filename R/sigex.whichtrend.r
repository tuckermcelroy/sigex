sigex.whichtrend <- function(mdl)
{

	##########################################################################
	#
	#	sigex.whichtrend
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
	#	Purpose: determines which component (index) corresponds to trend
	#	Background:
	#		A sigex model consists of process x = sum y, for 
	#		stochastic components y.  Each component process y_t
	#		is either stationary or is reduced to stationarity by
	#		application of a differencing polynomial delta(B), i.e.
	#			w_t = delta(B) y_t   is stationary.
	#		These differencing polynomials must be relatively prime.
	#		(This is an assumption of the code.)  If delta(1) = 0,
	#		then there is at least one unit root of that component,
	#		and hence that component corresponds to the stochastic trend.
	#		(It could have other non-stationary effects as well, but it
	#		at least contains some unit roots, and hence is designated
	#		as the trend component.)
	#	Inputs:
	#		mdl: the specified sigex model, a list object. 
	#			mdl[[1]] is mdlK, gives ranks of white noise covariance matrix
	#			mdl[[2]] is mdlType, a list giving t.s. model class, order, and bounds
	#			mdl[[3]] is mdlDiff, gives delta differencing polynomials
	#		      mdl[[4]] is list of regressors by individual series
	#	Outputs:
	#		trendcomp: index of the latent component that has a stochastic trend
	#
	####################################################################

	numcomp <- length(mdl[[3]])
	trendcomp <- NULL
	for(i in 1:numcomp)
	{
		if(sum(mdl[[3]][[i]])==0) trendcomp <- i
	}
	return(trendcomp)
}

