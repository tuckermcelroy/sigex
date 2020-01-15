#' converts a date in month-day-year format and finds the day of week
#'
#' @param date a given date in month-day-year format,
#'			a 3-element vector
#'
#' @return dayofweek: number 1 through 7, where 1 is Sunday, 2 is Monday, etc.
#'
#' @note calibrated according to Jan 1, 1600 being a Saturday
#'
#' @export
#'

day2week <- function(date)
{

	##########################################################################
	#
	#	day2week
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
	#	Purpose: converts a date in month-day-year format and finds the day of week
	#
	#	Inputs:
	#		date: a given date in month-day-year format,
	#			a 3-element vector
	#	Outputs:
	#		dayofweek: number 1 through 7, where 1 is Sunday, 2 is Monday, etc.
	#	Notes:
	#		calibrated according to Jan 1, 1600 being a Saturday
  # Require: date2day
	#
	#####################################################################

	month <- date[1]
	day <- date[2]
	year <- date[3]
	total.days <- 0
	if(year > 1600) {
	for(years in 1600:(year-1))
	{
		total.days <- total.days + date2day(12,31,years)
	} }
	total.days <- total.days + date2day(month,day,year) -1
	dayofweek <- total.days %% 7
	if(dayofweek==0) dayofweek <- 7
	return(dayofweek)
}

