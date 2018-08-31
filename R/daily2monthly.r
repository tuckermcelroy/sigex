#' converts a daily time series to a monthly time series
#'
#' @param daily.series the input daily time series
#' @param start.date beginning date in month-day-year format,
#'			a 3-element vector
#' @param dayofweek day of the week, 1 through 7
#'
#' @return monthly.series: each monthly value is the aggregate of
#'			all daily values in that month, of the specified
#'			day of the week
#' @export
#'


daily2monthly <- function(daily.series,start.date,dayofweek)
{

	##########################################################################
	#
	#	daily2monthly
	# 	    Copyright (C) 2018  Tucker McElroy
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
	#	Purpose: converts a daily time series to a monthly time series
	#
	#	Inputs:
	#		daily.series: the input daily time series
	#		start.date: beginning date in month-day-year format,
	#			a 3-element vector
	#		dayofweek: day of the week, 1 through 7
	#	Outputs:
	#		monthly.series: each monthly value is the aggregate of
	#			all daily values in that month, of the specified
	#			day of the week
	#	Requires:
	# 		day2date, date2day, day2week
	#
	#####################################################################

	T <- length(daily.series)
	my.dates <- t(apply(as.matrix(seq(0,T-1,1)),1,day2date,start.date))
	my.days <- rep(1,floor(T/7)+3) %x% t(seq(1,7))
	my.days <- matrix(t(my.days),ncol=1)
	my.days <- my.days[day2week(my.dates[1,]):(T+day2week(my.dates[1,])-1)]
 	dated.daily <- cbind(daily.series,my.days,my.dates)
	sub.series <- dated.daily[dated.daily[,2]==dayofweek,]
	monthly.series <- unique(sub.series[,c(3,5)])
	monthly.sums <- NULL
	for(i in 1:dim(monthly.series)[1])
	{
		indices <- (sub.series[,3]==monthly.series[i,1]) +
			(sub.series[,5]==monthly.series[i,2])
		new.sum <- sum(sub.series[indices==2,1])
		monthly.sums <- c(monthly.sums,new.sum)
	}
	monthly.series <- cbind(monthly.sums,monthly.series)

	return(monthly.series)
}


