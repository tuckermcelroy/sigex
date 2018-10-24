#' load data into a time series object
#'
#' @param data a T x N matrix, corresponding to N time series of length T
#' @param start.date date of first time obersvation; the
#'			 format is c(year,season)
#' @param period number of seasons per year
#' @param epithets vector of N character strings, giving a short name for
#'			 each series
#' @param plot boolean, whether to plot the series (max of N=10 allowed)
#'
#' @return data.ts: a T x N matrix ts object
#' @export
#'

sigex.load <- function(data,start.date,period,epithets,plot=FALSE)
{

	##########################################################################
	#
	#	sigex.load
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
	#	Purpose: load data into a time series object
	#
	#	Inputs:
	#		data: a T x N matrix, corresponding to N time series of length T
	#		start.date: date of first time obersvation; the
	#			 format is c(year,season)
	#		period: number of seasons per year
	#		epithets: vector of N character strings, giving a short name for
	#			 each series
	#		plot: boolean, whether to plot the series (max of N=10 allowed)
	#	Outputs:
	#		data.ts: a T x N matrix ts object
	#
	####################################################################

	data.ts <- ts(data,start=start.date,frequency=period,names=epithets)
 	if(plot) { plot(data.ts,xlab="Year") }

	return(data.ts)
}
