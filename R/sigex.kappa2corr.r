#' Transform N choose 2 kappa variables	into correlations, coherently
#' 
#' @param kappa a vector with N choose 2 real elements in [-1,1]
#' 
#' @return R.mat: an N x N correlation matrix 
#' @export
#' 

sigex.kappa2corr <- function(kappa)
{
  
  ##########################################################################
  #
  #	sigex.kappa2corr
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
  #	Purpose: transform N choose 2 kappa variables	into correlations, coherently.
  #		Based on bijection of Euclidean space and the set of correlation matrices.
  #
  #	Inputs:
  #		kappa: a vector with N choose 2 real elements in [-1,1]
  #	Outputs:
  #	  R.mat: an N x N correlation matrix 
  #
  #####################################################################
  
  N <- dim(kappa)[1]
  eta <- diag(N)
  for(i in 2:N)
  {
    for(j in 1:(i-1))
    {
      product <- 1
      if(j > 1) {
        for(l in 1:(j-1))
        {
          product <- product*sqrt(1 - kappa[i,l]^2)
        } }
      eta[i,j] <- kappa[i,j]*product
    }
    eta[i,i] <- product*sqrt(1 - kappa[i,i-1]^2)
  }
  R.mat <- eta %*% t(eta)
  
  return(R.mat)
}
