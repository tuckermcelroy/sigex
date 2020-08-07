#' obtain start and end dates for a weekly time series
#'
#' @param  first.day: for the weekly data, first.day indicates what day of
#'     the week corresponds to the first day in the week, 
#'     with 1 = Sunday, 2 = Monday, etc.
#'     e.g. for weekly data beginning with Sunday, set first.day = 1
#' @param  begin: two components with year and week index (between 1 and 53)
#' @param  T: number of weeks in the time series
#'
#' @return list with start.date and end.date, each in format c(month,day,year)
#' @export
#'

weekly2date <- function(first.day,begin,T)
{
  
  ##########################################################################
  #
  #	weekly2date
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
  #	Purpose: obtain start and end dates for a weekly time series
  #	Inputs:
  #   first.day: for the weekly data, first.day indicates what day of
  #     the week corresponds to the first day in the week, 
  #     with 1 = Sunday, 2 = Monday, etc.
  #     e.g. for weekly data beginning with Sunday, set first.day = 1
  #   begin: two components with year and week index (between 1 and 53)
  #   T: number of weeks in the time series
  #	Outputs:
  #		list with start.date and end.date, each in format c(month,day,year)
  # Require: date2day, day2week, day2date
  #
  #####################################################################
  
  start.year <- begin[1]
  start.week <- begin[2]
  day.lead <- day2week(c(1,1,start.year)) - first.day
  if(day.lead < 0) { day.lead <- day.lead + 7 }
  day.index <- 7*(start.week-1) - day.lead + 1
  year.index <- start.year
  if(day.index <= 0)
  {
    day.index <- date2day(12,31,start.year-1) + day.index
    year.index <- year.index-1
  }
  start.date <- day2date(day.index-1,c(1,1,year.index))
  end.date <- day2date(day.index-2 + 7*T,c(1,1,year.index))
  
  return(list(start.date,end.date))
}

