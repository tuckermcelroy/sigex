sigex.renderpd <- function(L.mat,D.mat,thresh)
{
	
	##########################################################################
	#
	#	sigex.renderpd
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
	#	Purpose: modifies a covariance matrix so that it is pd
	#	Background: a non-negative definite matrix Sigma has a
	#		Generalized Cholesky Decomposition (GCD) of the form
	#		Sigma = L %*% D %*% t(L),
	#		where L is unit lower triangular and D is diagonal with
	#		non-negative entries, referred to as the Schur complements
	#		of Sigma.  The number of nonzero Schur complements equals
	#		the rank of Sigma.   The condition numbers can be computed
	#		by dividing D by the diagonal of Sigma.
	#	Inputs:
	#		L.mat: rectangular, lower Cholesky factor
	#		D.mat: vector of diagonal entries in Cholesky decomposition
	#		thresh: a threshold for condition numbers.  The new 
	#			covariance matrix has a modified D matrix, such
	#			that the condition numbers (with the same L.mat)
	#			are bounded below by thresh.
	#	Note: if run on a reduced rank matrix, pad out
	#		D.mat with -Inf and corresponding columns
	#		of L.mat should be unit vector.
	#		Dimension >= 2, because first condition number is zero.
	#	Outputs:
	#		conds.new: new condition numbers
	#
	####################################################################

	N <- dim(L.mat)[1]
	D.new <- matrix(exp(D.mat[1]),1,1)
	for(i in 2:N)
	{
		val <- (matrix(L.mat[i,1:(i-1)],nrow=1) %*% D.new %*% matrix(L.mat[i,1:(i-1)],ncol=1))/(exp(-thresh) - 1)
		d.new <- max(exp(D.mat[i]),val)
		D.new <- diag(c(diag(D.new),d.new))
	}
	conds.new <- log(diag(D.new))

	return(conds.new)
}
