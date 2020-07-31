sigex.castextract <- function(data.ts,data.casts,mdl,castspan,param)
{
  
  ##########################################################################
  #
  #	sigex.castextract
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
  
  ################# Documentation ############################################
  #
  #	Purpose: formats cast estimates with two standard errors
  #	Background:	
  #		psi refers to a vector of real numbers containing all
  #		hyper-parameters (i.e., reals mapped bijectively to the parameter manifold)
  #	Inputs:
  #		data.ts: a T x N matrix ts object; any  values to be imputed
  #			must be encoded with NA in that entry.  The NA is for missing value,
  #     or an enforced imputation (e.g. extreme-value adjustment).
  #   data.casts: a list containing casts.x and casts.var 
  #	  	casts.x: N x H matrix of forecasts, midcasts, aftcasts, where H
  #		  	is the total number of time indices with missing values,
  #			  given by cardinality( leads setminus {1,2,...,T} )	
  #		  casts.var: NH x NH matrix of covariances of casting errors.
  #			  note that casts.var.array <- array(casts.var,c(N,H,N,H)) 
  #	  		corresponds to cast.var.array[,j,,k] equal to the 
  #		  	covariance between the jth and kth casting errors
  #		mdl: the specified sigex model, a list object
  #   castspan: an non-negative integer horizon giving number of fore- and aft-casts
  #		param: must have form specified by mdl
  #	Outputs:
  #		list object extract.casts with cast, upp, and low
  #		  cast: T x N matrix of the data or casts
  #		  upp: as cast, plus twice the standard error
  # 		low: as cast, minus twice the standard error
  # 
  ####################################################################
  
  x <- t(data.ts)
  N <- dim(x)[1]
  T <- dim(x)[2]
  tol <- 1e-10
  
  leads.fore <- NULL
  leads.aft <- NULL
  if(castspan > 0)
  {
    leads.fore <- seq(T+1,T+castspan)
    leads.aft <- seq(1-castspan,0)
  }
  if(length(leads.fore)>0) { x <- cbind(x,matrix(NA,nrow=N,ncol=length(leads.fore))) }
  if(length(leads.aft)>0) { x <- cbind(matrix(NA,nrow=N,ncol=length(leads.aft)),x) }
  TH <- dim(x)[2]
  
  times.na <- NULL
  for(k in 1:N)
  { times.na <- union(times.na,seq(1,TH)[is.na(x)[k,]]) }
  times.na <- sort(times.na)
  
  reg.trend <- ts(mdl[[4]][[1]][,1])
  d <- 0
  diff.flag <- TRUE
  while(diff.flag)
  {
    reg.trend <- ts(diff(reg.trend)[-1])
    if(sum(reg.trend^2)==0) { diff.flag <- FALSE } else { d <- d+1 }
  }
  reg.trend <- t(param[[4]]) %x% seq(1-castspan,T+castspan)^d
 
  x.casts <- data.casts[[1]] 
  se.casts <- matrix(sqrt(tol+diag(data.casts[[2]])),nrow=N)
  extract.casts <- list()
  extract.casts[[1]] <- t(x) - reg.trend
  extract.casts[[1]][times.na,] <- t(x.casts)
  extract.casts[[1]] <- extract.casts[[1]] + reg.trend
  extract.casts[[2]] <- t(x) - reg.trend
  extract.casts[[2]][times.na,] <- t(x.casts) + 2*t(se.casts)
  extract.casts[[2]] <- extract.casts[[2]] + reg.trend
  extract.casts[[3]] <- t(x) - reg.trend
  extract.casts[[3]][times.na,] <- t(x.casts) - 2*t(se.casts)
  extract.casts[[3]] <- extract.casts[[3]] + reg.trend
  
  return(extract.casts)
}
