#' compute the sum of two polynomials,
#		assuming both are symmetric power series
#'
#' @param a symmetric vector of polynomial coefficients,
#'			where a[1] is the zeroth coefficient,
#'			a[2] is the first coefficient, etc.
#' @param b symmetric vector of polynomial coefficients,
#'			where b[1] is the zeroth coefficient,
#'			b[2] is the first coefficient, etc.
#'
#' @return symmetric vector of polynomial coefficients for c(z) = a(z) + b(z),
#'			where c[1] is the zeroth coefficient, c[2] is the first coefficient, etc.
#' @export
#'
polysum <- function(a,b)
{

	##########################################################################
	#
	#	polysum
	# 	    Copyright (C) 2017  Tucker McElroy
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
	#	Purpose: compute the sum of two polynomials,
	#		assuming both are symmetric power series
	#	Inputs:
	#		a: symmetric vector of polynomial coefficients,
	#			where a[1] is the zeroth coefficient,
	#			a[2] is the first coefficient, etc.
	#		b: symmetric vector of polynomial coefficients,
	#			where b[1] is the zeroth coefficient,
	#			b[2] is the first coefficient, etc.
	#	Outputs:
	#		c: symmetric vector of polynomial coefficients for c(z) = a(z) + b(z),
	#			where c[1] is the zeroth coefficient, c[2] is the first coefficient, etc.
	#
	##########################################################################################

      n <- length(a)
      m <- length(b)
      if (m > n) out <- b + c(rep(0,(m-n)/2),a,rep(0,(m-n)/2))
      if (n > m) out <- a + c(rep(0,(n-m)/2),b,rep(0,(n-m)/2))
      if (n==m) out <- a + b
      return(out)
}
