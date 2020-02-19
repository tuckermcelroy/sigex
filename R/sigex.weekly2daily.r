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
  # Requires: day2week, date2day
  #
  ##############################################################
  
  start.week <- start(data.ts)[2]
  day.inyear <- 7*(start.week-1) + 1 + first.day
  
  #HERE
  
  start.day <- date2day(start.date[1],start.date[2],start.date[3])
  week.days <- c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday")
  days.index <- seq(first.day,first.day+6) %% 7
  days.index[days.index==0] <- 7
  ragged.fore <- day2week(start.date) - first.day
  if(ragged.fore > 0) { data.ts <- c(rep(NA,ragged.fore),data.ts) }
  ragged.aft <- length(data.ts) %% 7
  if(ragged.aft > 0) { data.ts <- c(data.ts,rep(NA,7-ragged.aft))}
  data.mat <- t(matrix(data.ts,nrow=7))
  data.ts <- ts(data.mat,start=c(start.date[3],ceiling(start.day/7)),frequency=52,
                names=week.days[days.index])
  
  return(data.ts)
}


