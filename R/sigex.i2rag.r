sigex.i2rag <- function(z)
{

  ##########################################################################
  #
  #	sigex.i2rag
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
  #	Purpose:  read multivariate time series and generate list of missing value indices
  #	Inputs:
  #		z: raw data as N x T matrix with missing values at various time points.
  #			Missing values are at any 1 <= t <= T, and occur for some of the N series
  #     (the ragged case), and are denoted with a 1i.  That is, 
  #			Im(z[,t]) = rep(1i,N) or subset thereof encodes missing values.
  #	Outputs:
  #		list containing leads.rag and ragged
  #		leads.rag: an integer sequence of indices that have at least one missing value
  #			These integers are a subset of {1,2,...,T}.  
  #   ragged: a list with number of items equal to length of leads.rag
  #     within the sample {1,2,...,T}.  (ragged=NULL if there are none.)
  #     Each element contains indices of vector components that are missing.
  #     Note: ragged[[j]] must be non-empty, where leads[j] is the time
  #      where at least one missing value (in-sample) occurs.
  #
  ####################################################################
 
  N <- dim(z)[1]
  T <- dim(z)[2]
  all.series <- seq(1,N)
  all.indices <- seq(1,T)
  full.indices <- all.indices[colSums(Im(z)==1)==0]
  cast.indices <- setdiff(all.indices,full.indices)
  ragged <- list()
  leads.rag <- NULL
  for(t in 1:length(cast.indices))
  {
    rag.series <- all.series[Im(z[,cast.indices[t]])==1,drop=FALSE]
    if(length(rag.series)<=N) 
    { 
      ragged[[length(ragged)+1]] <- rag.series 
      leads.rag <- c(leads.rag,cast.indices[t])
    }
  }

  return(list(leads.rag,ragged))
}

