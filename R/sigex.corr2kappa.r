#' Transform correlation matrix into N choose 2 kappa variables
#' 
#' @param corr.mat an N x N correlation matrix 
#' 
#' @return kappa: N choose 2 vector of [-1,1] numbers
#' @export
#' 

sigex.corr2kappa <- function(corr.mat)
{
  
  ##########################################################################
  #
  #	sigex.corr2kappa
  # 	    Copyright (C) 2024  Tucker McElroy
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
  #	Purpose: transform correlation matrix into N choose 2 kappa variables.
  #		Based on bijection of Euclidean space and the set of correlation matrices.
  #
  #	Inputs:
  #		corr.mat: an N x N correlation matrix 
  #	Outputs:
  #	  kappa: N choose 2 vector of [-1,1] numbers
  #
  #####################################################################
  
  N <- dim(corr.mat)[1]
  corr.gcd <- getGCD(corr.mat,N)
  L.mat <- corr.gcd[[1]]
  D.mat <- corr.gcd[[2]]
  rank.config <- seq(1,N)[D.mat > 10e-15]
  new.z <- diag(N)
  new.z[lower.tri(new.z)] <- NA
  L.mat <- L.mat[,rank.config]
  D.mat <- D.mat[rank.config]
  D.inv <- rep(1,N)
  D.inv[rank.config] <- 1/sqrt(D.mat)
  sub.z <- diag(D.inv) %*% L.mat %*% diag(sqrt(D.mat))
  new.z[,rank.config] <- sub.z
  w <- new.z
  w[,-rank.config] <- 0
  
  kappa <- diag(N)
  for(i in 2:N)
  {
    for(j in 1:(i-1))
    {
      summand <- w[i,i]^2
      for(l in j:(i-1))
      {
        summand <- summand + w[i,l]^2 
      }
      if(summand==0) { kappa[i,j] <- 0 } else {
        kappa[i,j] <- w[i,j]/sqrt(summand) }
    }
  }
  
  return(kappa)
}