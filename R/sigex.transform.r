#' applies aggregation, followed by transformations to the data
#'
#' @param data.ts a T x N matrix ts object,
#'			corresponding to N time series of length T
#' @param transform a character indicating an instantaneous
#'			transformation to be applied; current options are
#'			"none", "log", and "logistic"
#' @param aggregate a boolean, set to TRUE if all series are to
#'			be aggregated into a total
#'
#' @return data.ts: a T x N0 matrix ts object, where N0=1 if
#'			aggregate=TRUE, otherwise N0=N
#' @export
#'

sigex.transform <- function(data.ts,transform,aggregate=FALSE)
{
	##########################################################################
	#
	#	sigex.transform
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
	#	Purpose: applies aggregation, followed by transformations to the data
	#
	#	Inputs:
	#		data.ts: a T x N matrix ts object,
	#			corresponding to N time series of length T
	#		transform: a character indicating an instantaneous
	#			transformation to be applied; current options are
	#			"none", "log", and "logistic"
	#		aggregate: a boolean, set to TRUE if all series are to
	#			be aggregated into a total
	#	Outputs:
	#		data.ts: a T x N0 matrix ts object, where N0=1 if
	#			aggregate=TRUE, otherwise N0=N
	#
	####################################################################

	if(aggregate) {
		data <- as.matrix(rowSums(data.ts))
		new.names <- "Total"
	} else {
		data <- data.ts
		new.names <- colnames(data.ts)
	}
	if(transform=="log") data <- log(data)
	if(transform=="logistic") data <- log(data) - log(1-data)
	if(transform=="none") data <- data

	data.ts <- ts(data,start=start(data.ts),frequency=frequency(data.ts),
		names=new.names)

	return(data.ts)
}


