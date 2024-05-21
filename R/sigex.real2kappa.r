#'  Transform extended real variables into N choose 2 kappa variables
#'  
#'  @param z a vector with N choose 2 real elements	in [-Inf,Inf]
#'  
#'  @return kappa: kappa values in [-1,1], in lotri form
#'  @export
#'  

sigex.real2kappa <- function(z)
{
  
  ##########################################################################
  #
  #	sigex.real2kappa
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
  #	Purpose: transform extended real variables into N choose 2 kappa variables
  #
  #	Inputs:
  #		z: a vector with N choose 2 real elements	in [-Inf,Inf]
  #	Outputs:
  #	  kappa: kappa values in [-1,1], in lotri form
  #
  #####################################################################
  
  M <- length(z)
  N <- max(Re(polyroot(c(-2*M,-1,1))))
  for(j in 1:M) {
    if(z[j] == -Inf) { z[j] <- -1 } else { 
      z[j] <- (1 - exp(-1*z[j]))/(1 + exp(-1*z[j])) } }
  kappa <- diag(N)
  kappa[lower.tri(kappa)] <- z
  
  return(kappa)
}

