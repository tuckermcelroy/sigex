VARMAauto <- function(phi,theta,sigma,maxlag)
{

	##########################################################################
	#
	#	VARMAauto
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
	#	Purpose: computes autocovariances of VARMA
	#	Background: function computes autocovariances of VARMA (p,q) from lag zero
	#		to maxlag, with array inputs phi and theta.  VARMA equation:
	#	(1 - phi[1]B ... - phi[p]B^p) X_t = (1 + theta[1]B ...+ theta[q]B^q) WN_t
	#	Note: for absent VAR or VMA portions, pass in NULL
	#	Inputs:
	#		phi: array of dimension m x m x p of VAR coefficients, e.g.,
	#			phi <- array(cbind(phi1,phi2,...,phip),c(m,m,p))
	#		theta: array of dimension m x m x q of VMA coefficients, e.g.,
	#			theta <- array(cbind(theta1,theta2,...,thetaq),c(m,m,q))
	#		sigma: m x m covariance matrix of white noise
	#	Outputs:
	#		autocovariances at lags 0 through maxlag, as array of dimension m x m x (maxlag+1)
	#
	####################################################################

polymulMat <- function(amat,bmat)
{
	p <- dim(amat)[3]-1
	q <- dim(bmat)[3]-1
	N <- dim(amat)[2]

	r <- p+q
	bmat.pad <- array(0,c(N,N,r+1))
	for(i in 1:(q+1)) { bmat.pad[,,i] <- bmat[,,i] }
	cmat <- array(0,c(N,N,r+1))
	cmat[,,1] <- amat[,,1] %*% bmat.pad[,,1]
	for(j in 2:(r+1))
	{
		cmat[,,j] <- amat[,,1] %*% bmat.pad[,,j]
		for(k in 1:min(p,j-1))
		{ cmat[,,j] <- cmat[,,j] + amat[,,k+1] %*% bmat.pad[,,j-k] }
	}

	return(cmat)
}

Kcommut <- function(vect,m,n)
{
	return(matrix(t(matrix(vect,nrow=m,ncol=n)),ncol=1))
}

m <- dim(sigma)[2]
p <- 0
q <- 0
if (length(phi) > 0) p <- dim(phi)[3]
if (length(theta) > 0) q <- dim(theta)[3]
Kmat <- apply(diag(m^2),1,Kcommut,m,m)

if (q == 0) { gamMA <- array(sigma,c(m,m,1)) } else
{
	temp <- polymulMat(array(cbind(diag(m),matrix(theta,m,m*q)),c(m,m,q+1)),
		array(sigma,c(m,m,1)))
	flip <- aperm(theta,c(2,1,3))[,,q:1,drop=FALSE]
	gamMA <- polymulMat(temp,array(cbind(matrix(flip,m,m*q),diag(m)),c(m,m,q+1)))
}
gamMA <- gamMA[,,(q+1):(2*q+1),drop=FALSE]
gamMAvec <- matrix(gamMA,m^2*(q+1),1)

if (p > 0)
{
	Amat <- matrix(0,nrow=m^2*(p+1),ncol=m^2*(2*p+1))
	Amat <- array(Amat,c(m^2,p+1,m^2,2*p+1))
	Arow <- diag(m^2)
	for(i in 1:p)
	{
		Arow <- cbind(-1*diag(m) %x% phi[,,i],Arow)
	}
	for(i in 1:(p+1))
	{
		Amat[,i,,i:(i+p)] <- Arow
	}
	newA <- array(matrix(Amat[,1:(p+1),,1:p],m^2*(p+1),m^2*(p)),c(m^2,p+1,m^2,p))
	for(i in 1:(p+1))
	{
		for(j in 1:p)
		{
 			newA[,i,,j] <- newA[,i,,j] %*% Kmat
		}
	}
	Amat <- cbind(matrix(Amat[,,,p+1],m^2*(p+1),m^2),
			matrix(Amat[,,,(p+2):(2*p+1)],m^2*(p+1),m^2*(p)) +
			matrix(newA[,,,p:1],m^2*(p+1),m^2*(p)))

	Bmat <- matrix(0,nrow=m^2*(q+1),ncol=m^2*(p+q+1))
	Bmat <- array(Bmat,c(m^2,q+1,m^2,p+q+1))
	Brow <- diag(m^2)
	for(i in 1:p)
	{
		Brow <- cbind(Brow,-1*phi[,,i] %x% diag(m))
	}
	for(i in 1:(q+1))
	{
		Bmat[,i,,i:(i+p)] <- Brow
	}
	Bmat <- Bmat[,,,1:(q+1)]
	Bmat <- matrix(Bmat,m^2*(q+1),m^2*(q+1))
# modification suggested by A. Roy
#	Binv <- solve(Bmat)
#	gamMix <- Binv %*% gamMAvec
	gamMix <- solve(Bmat,gamMAvec)
	if (p <= q) gamMixTemp <- gamMix[1:((p+1)*m^2)] else
		gamMixTemp <- c(gamMix,rep(0,(p-q)*m^2))
	# gamARMA <- solve(Amat) %*% gamMixTemp
	gamARMA <- solve(Amat, gamMixTemp)
	gamMix <- array(matrix(gamMix,m,m*(q+1)),c(m,m,q+1))
	gamARMA <- array(matrix(gamARMA,m,m*(p+1)),c(m,m,p+1))
} else
{
	gamARMA <- array(gamMA[,,1],c(m,m,1))
	if (q == 0) { gamMix <- array(sigma,c(m,m,1)) } else
		gamMix <- gamMA[,,1:(q+1)]
}

if (maxlag <= p)
{
	gamARMA <- gamARMA[,,1:(maxlag+1)]
} else
{
	if (maxlag > q) gamMix <- array(cbind(matrix(gamMix,m,m*(q+1)),
		matrix(0,m,m*(maxlag-q))),c(m,m,(maxlag+1)))
	for(k in 1:(maxlag-p))
	{
		len <- dim(gamARMA)[3]
		acf <- gamMix[,,p+1+k]
		if (p > 0)
		{
			temp <- NULL
			for(i in 1:p)
			{
				temp <- rbind(temp,gamARMA[,,len-i+1])
			}
			acf <- acf + matrix(phi,m,m*p) %*% temp
		}
		gamARMA <- array(cbind(matrix(gamARMA,m,m*len),acf),c(m,m,len+1))
	}
}

return(gamARMA)
}
