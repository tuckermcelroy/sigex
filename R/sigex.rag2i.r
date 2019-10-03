sigex.rag2i <- function(leads.rag,ragged,z)
{
  
  ##########################################################################
  #
  #	sigex.rag2i
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
  #	Purpose:  read multivariate time series and insert 1i for indicated missing values
  #	Inputs:
  #		leads.rag: an integer sequence of indices that have at least one missing value
  #			These integers are a subset of {1,2,...,T}.  
  #   ragged: a list with number of items equal to length of leads.rag
  #     within the sample {1,2,...,T}.  (ragged=NULL if there are none.)
  #     Each element contains indices of vector components that are missing.
  #     Note: ragged[[j]] must be non-empty, where leads[j] is the time
  #      where at least one missing value (in-sample) occurs.
  #		z: raw data as N x T matrix with missing values at various time points.
  #			Missing values are at any 1 <= t <= T, and occur for some of the N series
  #     (the ragged case), and are denoted with a 1i.  That is, 
  #			Im(z[,t]) = rep(1i,N) or subset thereof encodes missing values.
  #	Outputs:
  #		z: like z, but with additional 1i inserted where indicated by ragged
  #
  ####################################################################
  
  N <- dim(z)[1]
  T <- dim(z)[2]
  if(length(leads.rag)>0) 
  { 
    for(t in 1:length(leads.rag))
    {
      raggeds <- ragged[[t]]
      z[raggeds,leads.rag[t]] <- rep(1i,length(raggeds)) 
    }
  }
  
  return(z)
}

  