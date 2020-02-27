sigex.weekly2daily <- function(data.ts,first.day)
{
  
  ##########################################################################
  #
  #	sigex.weekly2daily
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
  #	Purpose: de-embeds a weekly time series as a daily time series
  #	Inputs:
  #		data.ts: a T x 7 matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #   first.day is a number between 1 and 7 that indicates what day of the week
  #     should correspond to the first component of the 7-vector.
  #     So 1 = Sunday, 2 = Monday, 3 = Tuesday, 4 = Wednesday, 
  #       5 = Thursday, 6 = Friday, 7 = Sunday.
  #	outputs:
  #   data.ts: a T x 1 matrix ts object; any  values to be imputed
  #			are marked with NA in that entry.  The NA is for missing value.
  # Requires: day2week, date2day, sigex.load
  #
  ##############################################################
  
  start.year <- start(data.ts)[1]
  start.week <- start(data.ts)[2]

  day.lead <- day2week(c(1,1,start.year)) - first.day
  if(day.lead < 0) { day.lead <- day.lead + 7 }
  day.index <- 7*(start.week-1) - day.lead + 1
  year.index <- start.year
  if(day.index <= 0)
  {
    day.index <- date2day(12,31,start.year-1) + day.index
    year.index <- year.index-1
  }
  #start.date <- day2date(day.index,c(1,1,year.index))
  data.ts <- matrix(t(data.ts),ncol=1)
  data.ts <- sigex.load(data.ts,c(year.index,day.index),365,"",FALSE)
  
  return(data.ts)
}


