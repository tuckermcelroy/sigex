x11filters <- function(period,p.seas)
{
  
  ##########################################################################
  #
  #	x11filters
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
  #	Purpose: generates x11 trend, seasonal, and seasonal adjustment filters.
  #	Inputs:	
  #   period: possibly non-integer period of seasonality
  #   p.seas: integer number of seasonal moving averages
  #	Outputs: 
  #	  list of arrays for trend.filter, seas.filter, and sa.filter
  # Requires: ubgenerator, polymult
  #
  ##########################################################################################
  
  half.len <- floor(period/2)
  p.seas <- 1
  trend.filter <- ubgenerator(period,NULL,1000)
  trend.filter <- trend.filter/sum(trend.filter)
  detrend.filter <- c(rep(0,half.len),1,rep(0,half.len)) - trend.filter
  seas.filter <- 0
  for(j in 1:p.seas)
  {
    periodj <- j*period
    half.lenj <- floor(periodj/2)
    seas.filterj <- ubgenerator(periodj,half.lenj-1,1000)
    seas.filterj <- polymult(seas.filterj,c(1,0,-1))
    seas.filter <- c(seas.filter,rep(0,length(seas.filterj)-length(seas.filter)))
    seas.filter <- seas.filter + seas.filterj
  }
  seas.filter <- c(rep(0,length(seas.filter)-1),seas.filter)
  seas.filter <- seas.filter + rev(seas.filter)
  factor <- 2*p.seas + 1
  seas.filter <- seas.filter/factor
  seas.filter <- c(rep(0,(length(seas.filter)-1)/2),1,rep(0,(length(seas.filter)-1)/2)) - seas.filter
  sa.filter <- polymult(detrend.filter,seas.filter)
  shift <- (length(sa.filter)-1)/2
  sa.filter <- c(1,rep(0,shift)) - rev(sa.filter[1:(shift+1)])
  sa.filter <- c(rev(sa.filter),sa.filter[-1])
  trend.filter <- array(trend.filter,c(1,1,length(trend.filter)))
  seas.filter <- array(seas.filter,c(1,1,length(seas.filter)))
  sa.filter <- array(sa.filter,c(1,1,length(sa.filter)))
  
  return(list(trend.filter,seas.filter,sa.filter))
}