#' plot AR spectrum of input time series
#'
#' @param data.ts a T x N matrix ts object,
#'			corresponding to N time series of length T
#' @param diff boolean, if TRUE difference the series
#'			and plot Growth Rate, else in Levels
#' @param subseries index between 1 and N, indicating which series
#'			to examine
#' @param period number of observations per cycle (e.g. a year or week)
#'
#' @return NA
#' @export
#'

sigex.specar <- function(data.ts,diff=FALSE,subseries,period)
{

	##########################################################################
	#
	#	sigex.specar
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
	#	Purpose: plot AR spectrum of input time series
	#
	#	Inputs:
	#		data.ts: a T x N matrix ts object,
	#			corresponding to N time series of length T
	#		diff: boolean, if TRUE difference the series
	#			and plot Growth Rate, else in Levels
	#		subseries: index between 1 and N, indicating which series
	#			to examine
	#		period: number of observations per cycle (e.g. a year or week)
	#	Outputs:
	#		NA
	#
	####################################################################

	data.ts <- ts(data.ts,frequency=period)
	freqs <- floor(period/2)
	if(diff) {
		spec.ar(diff(data.ts)[,subseries],
			main= paste(colnames(data.ts)[subseries],"Growth Rate"))
		abline(v=seq(1,freqs),col=2)
	} else {
		spec.ar(data.ts[,subseries],
			main= paste(colnames(data.ts)[subseries],"Levels"))
		abline(v=seq(1,freqs),col=2)
	}

}
