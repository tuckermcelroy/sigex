sigex.daily2weekly <- function(data.ts,first.day,start.date)
{
  
  ##########################################################################
  #
  #	sigex.daily2weekly
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
  #	Purpose: embeds a daily time series as a weekly time series
  #	Inputs:
  #		data.ts: a T x 1 matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #   first.day is a number between 1 and 7 that indicates what day of the week
  #     should correspond to the first component of the 7-vector.
  #     So 1 = Sunday, 2 = Monday, 3 = Tuesday, 4 = Wednesday, 
  #       5 = Thursday, 6 = Friday, 7 = Sunday.
  #		start.date: a given starting date in month-day-year format,
  #			a 3-element vector
  #	outputs:
  #   data.ts: a T x 7 matrix ts object; any  values to be imputed
  #			are marked with NA in that entry.  The NA is for missing value.
  # Requires: day2week, date2day
  #
  ##############################################################
  
  start.day <- date2day(start.date[1],start.date[2],start.date[3])
  week.days <- c("Sunday","Monday","Tuesday","Wednesday","Thursday","Friday","Saturday")
  days.index <- seq(first.day,first.day+6) %% 7
  days.index[days.index==0] <- 7
  ragged.fore <- day2week(start.date) - first.day
  if(ragged.fore < 0) { ragged.fore <- ragged.fore + 7 }
  if(ragged.fore > 0) { data.ts <- c(rep(NA,ragged.fore),data.ts) }
  ragged.aft <- length(data.ts) %% 7
  if(ragged.aft > 0) { data.ts <- c(data.ts,rep(NA,7-ragged.aft))}
  week.index <- ceiling((start.day - ragged.fore -1)/7) + 1
  data.mat <- t(matrix(data.ts,nrow=7))
  data.ts <- ts(data.mat,start=c(start.date[3],week.index),frequency=53,
                names=week.days[days.index])
  
  return(data.ts)
}

  
  