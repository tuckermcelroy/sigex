sigex.param2gcd <- function(L.psi,N,vrank)
{

	##########################################################################
	#
	#	sigex.param2gcd
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
	#	Purpose: utility that takes a real vector and inserts as
	#		entries of a unit lower triangular matrix
	#	Background: a non-negative definite matrix Sigma has a
	#		Generalized Cholesky Decomposition (GCD) of the form
	#		Sigma = L %*% D %*% t(L),
	#		where L is unit lower triangular and D is diagonal with
	#		non-negative entries, referred to as the Schur complements
	#		of Sigma.  The number of nonzero Schur complements equals
	#		the rank of Sigma.
	#	Inputs:
	#		L.psi: a vector of reals of length <= N(N-1)/2
	#		N: dimension of the time series
	#		vrank: vector of integers between 1 and N, corresponding
	#			to indices of non-zero Schur complements in the GCD
	#	Outputs:
	#		L.mat: unit lower	triangular matrix of dimension N x length(vrank) 
	#		with entries given by L.psi
	#
	####################################################################

	L.mat <- 1
	if(N > 1) {
	fill <- NULL
	ind <- 0
	for(j in vrank)
	{
		pad <- NULL
		if(j < N) pad <- L.psi[(ind+1):(ind + N-j)]
		new <- c(rep(0,j-1),1,pad)
		fill <- c(fill,new)
		ind <- ind + N - j
	}
	L.mat <- matrix(fill,nrow=N,ncol=length(vrank)) 
	}

	return(L.mat)
}

