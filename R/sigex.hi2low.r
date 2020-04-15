sigex.hi2low <- function(filter.hi,hi.freq,low.freq,shift.hi)
{
  
  ##########################################################################
  #
  #	sigex.hi2low
  # 	    Copyright (C) 2019  Tucker McElroy
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
  #	Purpose: embeds a given high frequency filter as a low frequency filter
  #     for a vector embedding
  #	Background:	a scalar high frequency filter acts on a high frequency time series
  #     by convolution.  If we embed as a low frequency time series, the
  #     corresponding filter is now a matrix filter.  
  #	Inputs:
  #  filter.hi is given scalar high frequency filter of length L
  #  hi.freq is sampling frequency of high frequency time series
  #  low.freq is sampling frequency of low frequency time series
  #  shift.hi gives the integer offset for the high frequency filter:
  #     filter coefficients have indices -shift.hi,...,0,...,L-1-shift.hi
  #     set shift.hi = 0 for a causal filter
  #	outputs:
  #   filter.low is array s x s x M,
  #     where s = hi.freq/low.freq (must be an integer) 
  #     and M is length, which depends on s, L, and shift.hi
  #   shift.low is integer offset for the low frequency filter
  ##############################################################
  
  s.embed <- hi.freq/low.freq
  len.hi <- length(filter.hi)
  if(shift.hi %% s.embed == 0) { left.pad <- NULL } else {
    left.pad <- rep(0,s.embed - (shift.hi %% s.embed)) }
  new.len <- length(c(left.pad,filter.hi))
  if(new.len %% s.embed == 0) { right.pad <- NULL } else {
    right.pad <- rep(0,s.embed - (new.len %% s.embed)) }
  right.pad <- c(right.pad,rep(0,s.embed))
  next.column <- c(left.pad,filter.hi,right.pad)
  filter.embed <- next.column
  len.low <- length(next.column)/s.embed
  for(k in 1:s.embed)
  {
    next.column <- c(0,rev(rev(next.column)[-1]))
    filter.embed <- cbind(filter.embed,next.column)
  }
  filter.low <- array(filter.embed,c(s.embed,len.low,s.embed))
  filter.low <- aperm(filter.low,c(1,3,2))
  shift.low <- (length(left.pad) + shift.hi)/s.embed
  
  return(list(filter.low,shift.low))
}


