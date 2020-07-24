#' computes the generalized Cholesky decomposition
#'
#' @param Sigma symmetric, non-negative definite matrix
#' @param Rank presumed rank of Sigmma, less than or equal
#'	  		to the dimension of the matrix
#'
#' @return list of length 2 giving L.mat and D.mat.
#'    L.mat: rectangular, lower Cholesky factor
#'		D.mat: vector of diagonal entries in Cholesky decomposition
#' @export
#'

getGCD <- function(Sigma,Rank)
{

	##########################################################################
	#
	#	getGCD
	# 	    Copyright (C) 2018  Tucker McElroy
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
	#	Purpose: computes the generalized Cholesky decomposition
	#
	#	Inputs:
	#		Sigma: symmetric, non-negative definite matrix
	#		Rank: presumed rank of Sigma, less than or equal
	#	  		to the dimension of the matrix
	#	Outputs:
	#		L.mat: rectangular, lower Cholesky factor
	#		D.mat: vector of diagonal entries in Cholesky decomposition
	#	Notes: Sigma = L D L', where D is diagonal of dimension equal to
	#		the rank, with positive entries, and L is unit lower triangular
	#		with number of columns equal to rank
	#
	#####################################################################

	N <- dim(Sigma)[1]
	L.mat <- matrix(1,1,1)
	L.mat.inv <- L.mat
	D.mat <- Sigma[1,1]
	if(N > 1) {
	for(j in 2:N)
	{

		D.inv <- 1/D.mat
		D.inv[D.mat==0] <- 0
		new.sigma <- Sigma[j,1:(j-1)]
		if(j==2) { new.sigma <- as.matrix(new.sigma); L.mat <- as.matrix(L.mat) }
		new.l <- new.sigma %*% t(L.mat.inv)*D.inv
		new.l.tilde <- new.l %*% L.mat.inv
		L.mat <- cbind(L.mat,rep(0,(j-1)))
		L.mat <- rbind(L.mat,c(new.l,1))
		L.mat.inv <- cbind(L.mat.inv,rep(0,j-1))
		L.mat.inv <- rbind(L.mat.inv,c(-1*new.l.tilde,1))
		if(j==2) new.d <- Sigma[2,2] - new.l^2*D.mat
		if(j > 2) new.d <- Sigma[j,j] - new.l %*% diag(D.mat) %*% t(new.l)
		if(new.d <= 0) { new.d <- 0 }
		D.mat <- c(D.mat,new.d)
	} }

	rank.index <- rank(D.mat,ties.method="first")
	dims <- seq(1,N)[rank.index > (N-Rank)]

	L.mat <- matrix(L.mat[,dims],nrow=N,ncol=length(dims))
	D.mat <- D.mat[dims]

	return(list(L.mat,D.mat))
}



