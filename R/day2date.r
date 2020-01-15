#' converts a day index to a date in month-day-year format
#'
#' @param day integer, the number of days past the given start.date
#'			for which the date is desired; negative values reckon
#'			into the past
#' @param start.date a given starting date in month-day-year format,
#'			a 3-element vector
#'
#' @return end.date: the ending date  in month-day-year format,
#'			a 3-element vector
#' @export
#'

day2date <- function(day,start.date)
{

	##########################################################################
	#
	#	day2date
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
	#	Purpose: converts a day index to a date in month-day-year format
	#
	#	Inputs:
	#		day: integer, the number of days past the given start.date
	#			for which the date is desired; negative values reckon
	#			into the past
	#		start.date: a given starting date in month-day-year format,
	#			a 3-element vector
	#	Outputs:
	#		end.date: the ending date  in month-day-year format,
	#			a 3-element vector
	# Requires:  date2day
  #
	#####################################################################

	year.counter <- start.date[3]

	day.index <- day + date2day(start.date[1],start.date[2],start.date[3])
	day.counting <- day.index

	if(day.index > 0) {
	while(day.index > 0)
	{
		leap.flag <- 0
		if(((year.counter %% 4 == 0) && (year.counter %% 100 != 0)) ||
			(year.counter %% 400 == 0)){ leap.flag <- 1 }
		day.counting <- day.index
		day.index <- day.index - 365
		if(leap.flag) day.index <- day.index -1
		year.counter <- year.counter + 1
	}
	year.counter <- year.counter - 1
	day.index <- day.counting  }

	if(day.index <= 0) {
	year.counter <- year.counter - 1
 	while(day.index <= 0)
	{
		leap.flag <- 0
		if(((year.counter %% 4 == 0) && (year.counter %% 100 != 0)) ||
			(year.counter %% 400 == 0)){ leap.flag <- 1 }
		day.counting <- day.index
		day.index <- day.index + 365
		if(leap.flag) day.index <- day.index +1
		year.counter <- year.counter - 1
	}
	year.counter <- year.counter + 1  }

	if((1 <= day.index) && (day.index <= 31)) month.counter <- 1
	if((32 <= day.index) && (day.index <= (59+leap.flag))) month.counter <- 2
	if(((60+leap.flag) <= day.index) && (day.index <= (90+leap.flag))) month.counter <- 3
	if(((91+leap.flag) <= day.index) && (day.index <= (120+leap.flag))) month.counter <- 4
	if(((121+leap.flag) <= day.index) && (day.index <= (151+leap.flag))) month.counter <- 5
	if(((152+leap.flag) <= day.index) && (day.index <= (181+leap.flag))) month.counter <- 6
	if(((182+leap.flag) <= day.index) && (day.index <= (212+leap.flag))) month.counter <- 7
	if(((213+leap.flag) <= day.index) && (day.index <= (243+leap.flag))) month.counter <- 8
	if(((244+leap.flag) <= day.index) && (day.index <= (273+leap.flag))) month.counter <- 9
	if(((274+leap.flag) <= day.index) && (day.index <= (304+leap.flag))) month.counter <- 10
	if(((305+leap.flag) <= day.index) && (day.index <= (334+leap.flag))) month.counter <- 11
	if(((335+leap.flag) <= day.index) && (day.index <= (365+leap.flag))) month.counter <- 12
	day.counter <- day.index - date2day(month.counter,1,year.counter) + 1

	end.date <- c(month.counter,day.counter,year.counter)
	return(end.date)
}
