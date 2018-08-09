specFactmvar <- function(xAcf)
{

	##########################################################################
	#
	#	specFactmvar
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
	#	Purpose: compute the spectral factorization of a given multivariate p.d. sequence
	#	Background: any sequence of matrices with the property that A[j] = t(A[-j])
	#		has a Hermitian Fourier transform.  If in addition the eigenvalues
	#		of this FT, for each frequency, are all non-negative, then 
	#		the sequence is positive definite (p.d.),
	#		and can be factored as the product of a causal matrix power series
	#		with its conjugate transpose.   This power series is  called the spectral factorization.  
	#		In the case that the p.d. sequence has finitely many terms, 
	#		it corresponds to the acf of a VMA process and the VMA polynomial
	#		is given by the spectral factorization, once it is normalized such that
	#		the leading coefficient is an identity matrix (standard form for a VMA)
	#	Inputs:	
	#		xAcf: an array containing the given p.d. sequence.  It has
	#			dimension N x (q+1) x N, where xAcf[,1,] corresponds
	#			to the zeroth coefficient matrix, xAcf[,2,] is the first
	#			coefficient matrix, etc.
	#	Outputs: 
	#		theta: the spectral factorization as a list object
	#			theta[[1]] is an array of dimension N x N x (q+1),
	#			consisting of the VMA coefficients from index q through index zero
	#			(Note that theta[[1]][,,q+1] should be an identity matrix.)
	#			theta[[2]] is the white noise covariance matrix
	#			theta[[3]] tracks the number of iterations used in this algorithm,
	#				which is based on Bauer's method, as described in McElroy (2018, JTSA)
	#
 	##########################################################################################

	N <- dim(xAcf)[1]
	q <- dim(xAcf)[2] - 1
	Zzero <- matrix(0,nrow=N,ncol=N)
	eps <- 1
	thresh <- 10^(-15)
	Linvt <- diag(N)
	Dinv <- solve(xAcf[,1,])
	sqrtLam <- svd(xAcf[,1,])
	sqrtLam <- sqrtLam$u %*% diag(sqrt(sqrtLam$d),nrow=N) %*% t(sqrtLam$u)
	Gseq <- solve(sqrtLam)
	oldLam <- xAcf[,1,]
	Dsqrt <- sqrtLam
	gamSeq <- NULL
	count <- 1
	maxcount <- 500
	while((eps > thresh)&&(count<maxcount))
	{
		if(count <= q) nextgam <- xAcf[,count+1,] else nextgam <- Zzero
		gamSeq <- cbind(nextgam,gamSeq)
		Bseq <- gamSeq %*% Gseq
		Lam <- xAcf[,1,] - Bseq %*% t(Bseq)
		sqrtLam <- svd(Lam)
		sqrtLam <- sqrtLam$u %*% diag(sqrt(sqrtLam$d),nrow=N) %*% t(sqrtLam$u)
		Gseq <- rbind(cbind(Gseq,-1*Gseq%*%t(Bseq)%*%solve(sqrtLam)),
			cbind(t(rep(1,count) %x% Zzero),solve(sqrtLam)))
		count <- count+1
		if(count > q) eps <- sum((oldLam - Lam)^2)
		oldLam <- Lam
#		print(c(count,eps))
#		print(Lam)
	}
	Thetas <- cbind(Bseq,sqrtLam) %*% (diag(count) %x% solve(sqrtLam))
	arrayThetas <- array(Thetas,dim=c(N,N,count))
	out <- list(arrayThetas[,,(count-q):count],Lam,count)
	return(out)
}
