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

datefinder <- function(dayofweek,week,month,year)
{
  
  ##########################################################################
  #
  #	datefinder
  # 	    Copyright (C) 2020  Tucker McElroy
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
  #	Purpose: finds the date corresponding to a given holiday 
  #     described as "the [dayofweek] occuring in [week] of [month] of [year]"
  #
  #	Inputs:
  #		dayofweek: number 1 through 7, where 1 is Sunday, 2 is Monday, etc.
  #   week: number 1 through 5, for the week within the month
  #   month: number 1 through 12, for the month within the year
  #   year: integer year
  #	Outputs:
  #		date: a given date in month-day-year format, a 3-element vector
  # Note: if week = 5, the program gives the last dayofweek that occurs
  #   in that month, which could occur on week 4 or week 5!
  # Require: date2day, day2week, day2date
  #
  #####################################################################
  
  daycount <- date2day(month,1,year) -1
  if(month < 12) {
    lom <- date2day(month+1,1,year) - date2day(month,1,year) } else {
    lom <- date2day(1,1,year+1) - date2day(month,1,year) }
  daydiff <- dayofweek - day2week(c(month,1,year)) 
  if(daydiff < 0) { daydiff <- daydiff + 7 }
  inmonth <- daydiff + 7*(week-1)
  if(inmonth >= lom) { inmonth <- inmonth - 7 }
  daycount <- daycount + inmonth
  date <- day2date(daycount,c(1,1,year))
  return(date)
}

