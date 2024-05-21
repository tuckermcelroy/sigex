#' Transform N choose 2 kappa variables	into extended real variables
#' 
#' @param kappa a vector with N choose 2 real elements	in [-1,1]
#' 
#' @return z: reals in [-Inf,Inf]
#' @export
#' 

sigex.kappa2real <- function(kappa)
{
  
  ##########################################################################
  #
  #	sigex.kappa2real
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
  #	Purpose: transform N choose 2 kappa variables	into extended real variables.
  #
  #	Inputs:
  #		kappa: a vector with N choose 2 real elements	in [-1,1]
  #	Outputs:
  #	  z: reals in [-Inf,Inf]
  #
  #####################################################################
  
  z <- kappa[lower.tri(kappa)]
  z <- (log(1+z) - log(1-z))
  
  return(z)
}

