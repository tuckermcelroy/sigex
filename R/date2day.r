#' converts month-day-year date, and returns day index
#'
#' @param month 1 through 12, corresponding to the date's month
#' @param day 1 through 31, corresponding to the date's day
#' @param year 0 through ???, corresponding to the date's year (A.D.)
#'
#' @return day.count: day index within a year of the given date, from 1 to 366
#' @export
#'

date2day <- function(month,day,year)
{

	##########################################################################
	#
	#	date2day
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
	#	Purpose: converts month-day-year date, and returns day index
	#
	#	Inputs:
	#		month: 1 through 12, corresponding to the date's month
	#		day: 1 through 31, corresponding to the date's day
	#		year: 0 through ???, corresponding to the date's year (A.D.)
	#	Outputs:
	#		day.count: day index within a year of the given date, from 1 to 366
	#
	#####################################################################

	leap.flag <- 0
	if(((year %% 4 == 0) && (year %% 100 != 0)) || (year %% 400 == 0)){
		leap.flag <- 1}

	day.count <- 0
	if(month==1) { day.count <- day } else {
	for(month.count in 1:(month-1))
	{
		if(month.count==1) day.add <- 31
		if(month.count==2){
			day.add <- 28
			if(leap.flag) day.add <- 29
		}
		if(month.count==3) day.add <- 31
		if(month.count==4) day.add <- 30
		if(month.count==5) day.add <- 31
		if(month.count==6) day.add <- 30
		if(month.count==7) day.add <- 31
		if(month.count==8) day.add <- 31
		if(month.count==9) day.add <- 30
		if(month.count==10) day.add <- 31
		if(month.count==11) day.add <- 30
		day.count <- day.count + day.add
	}
	day.count <- day.count + day
	}

	return(day.count)
}
