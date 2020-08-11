ubgenerator <- function(period,trunc.len,m)
{
  
  ##########################################################################
  #
  #	ubgenerator
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
  #	Purpose: compute the product of unit root differencing operators.
  # Background: the goal is to annihilate periodic components of
  #   frequency 2*pi*k/period for k = 1, 2, ..., [period/2].
  #   Compute product_k  (1 - 2 cos(2*pi*k/period) z + z^2) by using
  #   cepstral method for higher accuracy.
  #	Inputs:	
  #   period: possibly non-integer period of sinusoids
  #   trunc.len: gives the maximal index of k, and is set to [period/2]
  #     when trunc.len = NULL
  #   m: large integer giving number of cepstral coefficients
  #	Outputs: 
  #	  wolds: coefficients of product polynomial.
  #
  ##########################################################################################

  ceps2wold <- function(ceps,q)
  {
    m <- length(ceps)
    if(q > m) {	ceps <- c(ceps,rep(0,q-m)) }
    wold <- 1
    wolds <- wold
    for(j in 1:q)
    {
      wold <- sum(seq(1,j)*ceps[1:j]*wolds[j:1])/j
      wolds <- c(wolds,wold)
    }
    return(wolds)
  }
  
  half.len <- floor(period/2)
  if(length(trunc.len)==0) { trunc.len <- half.len }
  ceps <- rep(0,m)
  
  for(ell in 1:m)
  {
    ceps[ell] <- -2*sum(cos(2*pi*ell*seq(1,trunc.len)/period))/ell
  }
  wolds <- ceps2wold(ceps,2*trunc.len)
  
  return(wolds)
}
